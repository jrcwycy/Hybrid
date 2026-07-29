"""
Standalone worker: runs ONE btrack config (motion-model-only, track() no
optimize()) in its own process and prints a single JSON line of results.

Run directly for a quick check:
    python run_single_ablation.py --npy_path well.npy --config_path cfg.json \
        --param MotionModel.accuracy --value 2.0 --n_frames 15 \
        --vol_x 1600 --vol_y 1200

Designed to be called via subprocess from ablation_driver.py so each config
is fully memory-isolated -- a segfault/OOM here just gives a nonzero
returncode in the parent, it doesn't kill your kernel.
"""

import sys
import json
import copy
import argparse

import numpy as np
import pandas as pd
import btrack


def set_nested(cfg: dict, dotted_key: str, value):
    """e.g. 'MotionModel.accuracy' -> cfg['TrackerConfig']['MotionModel']['accuracy'] = value"""
    parts = dotted_key.split(".")
    node = cfg["TrackerConfig"]
    for p in parts[:-1]:
        node = node[p]
    node[parts[-1]] = value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--objects_pkl_path", default=None,
                     help="pickle of an already-built list of PyTrackObject -- "
                          "FASTEST path, skips regionprops AND dataframe parsing entirely.")
    ap.add_argument("--npy_path", default=None,
                     help="raw segmentation array -- SLOW path, reruns regionprops every call")
    ap.add_argument("--df_path", default=None,
                     help="precomputed segment dataframe (parquet/csv) -- FAST path, "
                          "skips regionprops entirely. Use this if you already have "
                          "ID/x/y/frame(/features) extracted.")
    ap.add_argument("--well_id", default=None,
                     help="required if --df_path is used and df contains multiple wells")
    ap.add_argument("--well_col", default="well")
    ap.add_argument("--id_col", default="ID")
    ap.add_argument("--frame_col", default="frame")
    ap.add_argument("--x_col", default="x")
    ap.add_argument("--y_col", default="y")
    ap.add_argument("--config_path", required=True)
    ap.add_argument("--param", default=None, help="dotted key, e.g. MotionModel.accuracy")
    ap.add_argument("--value", default=None, help="value to set (parsed as float if possible)")
    ap.add_argument("--n_frames", type=int, default=None,
                     help="subset to first N frames for faster/lighter iteration")
    ap.add_argument("--vol_x", type=int, required=True)
    ap.add_argument("--vol_y", type=int, required=True)
    ap.add_argument("--features", default="", help="comma-separated feature names, optional")
    ap.add_argument("--tracking_updates", default="MOTION",
                     help="comma-separated, e.g. MOTION or MOTION,VISUAL")
    ap.add_argument("--crop", default=None,
                     help="x0,x1,y0,y1 -- if set, keep only objects with x,y inside this "
                          "box (tests whether candidate-pool size per frame matters)")
    ap.add_argument("--update_method", default="EXACT", choices=["EXACT", "APPROXIMATE"])
    ap.add_argument("--search_radius", type=float, default=None,
                     help="only used by APPROXIMATE update method")
    ap.add_argument("--optimize", action="store_true",
                     help="also run tracker.optimize() (hypothesis model / global stitching)")
    args = ap.parse_args()

    features = tuple(f for f in args.features.split(",") if f)

    if args.objects_pkl_path is not None:
        # FASTEST path: object DATA already extracted (plain dicts), just
        # unpickle and rebuild PyTrackObjects. Note: PyTrackObject itself
        # can't be pickled (ctypes pointers), so this expects a pickle of
        # plain dicts -- see objects_to_dicts() / save_objects_pkl() helpers.
        import pickle
        with open(args.objects_pkl_path, "rb") as f:
            obj_dicts = pickle.load(f)
        if args.n_frames is not None:
            obj_dicts = [d for d in obj_dicts if d["t"] < args.n_frames]
        objects = [btrack.btypes.PyTrackObject.from_dict(d) for d in obj_dicts]

    elif args.df_path is not None:
        # FAST path: build objects directly from the precomputed dataframe,
        # no regionprops recomputation.
        if args.df_path.endswith(".parquet"):
            df = pd.read_parquet(args.df_path)
        else:
            df = pd.read_csv(args.df_path)

        if args.well_id is not None:
            df = df[df[args.well_col] == args.well_id]

        if args.n_frames is not None:
            keep_frames = sorted(df[args.frame_col].unique())[: args.n_frames]
            df = df[df[args.frame_col].isin(keep_frames)]

        objects = []
        for _, row in df.iterrows():
            obj_dict = {
                "ID": int(row[args.id_col]),
                "t": int(row[args.frame_col]),
                "x": float(row[args.x_col]),
                "y": float(row[args.y_col]),
                "z": 0.0,
                "label": 5,
            }
            for feat in features:
                obj_dict[feat] = float(row[feat])
            objects.append(btrack.btypes.PyTrackObject.from_dict(obj_dict))

    elif args.npy_path is not None:
        # SLOW path: recompute from raw segmentation every call.
        seg = np.load(args.npy_path)
        if args.n_frames is not None:
            seg = seg[: args.n_frames]
        objects = btrack.utils.segmentation_to_objects(
            seg,
            properties=features if features else None,
            num_workers=4,
        )
    else:
        raise ValueError("Must supply either --df_path (fast) or --npy_path (slow)")

    if args.crop is not None:
        x0, x1, y0, y1 = (float(v) for v in args.crop.split(","))
        objects = [o for o in objects if x0 <= o.x <= x1 and y0 <= o.y <= y1]
        vol = ((x0, x1), (y0, y1))
    else:
        vol = ((0, args.vol_x), (0, args.vol_y))

    with open(args.config_path) as f:
        cfg = json.load(f)

    if args.param is not None:
        try:
            value = json.loads(args.value)
        except (json.JSONDecodeError, TypeError):
            value = args.value
        set_nested(cfg, args.param, value)

    tmp_cfg_path = "/tmp/_worker_cfg.json"
    with open(tmp_cfg_path, "w") as f:
        json.dump(cfg, f)

    tracking_updates = args.tracking_updates.split(",")

    with btrack.BayesianTracker() as tracker:
        tracker.configure(tmp_cfg_path)
        tracker.tracking_updates = tracking_updates
        if features:
            tracker.features = list(features)
        tracker.append(objects)
        tracker.volume = vol
        tracker.update_method = getattr(btrack.constants.BayesianUpdates, args.update_method)
        if args.search_radius is not None:
            tracker.max_search_radius = args.search_radius
        print(f"DEBUG: n_objects={len(objects)} update_method={tracker.update_method} "
              f"max_search_radius={getattr(tracker, 'max_search_radius', None)}",
              file=sys.stderr)
        tracker.track()
        if args.optimize:
            tracker.optimize()
        lengths = np.array([len(t) for t in tracker.tracks])

    result = {
        "param": args.param,
        "value": args.value,
        "n_tracks": int(len(lengths)),
        "mean_len": float(lengths.mean()) if len(lengths) else 0.0,
        "median_len": float(np.median(lengths)) if len(lengths) else 0.0,
        "frac_len1": float((lengths == 1).mean()) if len(lengths) else 1.0,
        "frac_len_lt5": float((lengths < 5).mean()) if len(lengths) else 1.0,
        "frac_len_ge20": float((lengths >= 20).mean()) if len(lengths) else 0.0,
        "max_len": int(lengths.max()) if len(lengths) else 0,
    }
    # single JSON line on its own -- driver parses the LAST line of stdout
    print("RESULT_JSON:" + json.dumps(result))


if __name__ == "__main__":
    main()