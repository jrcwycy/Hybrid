"""
Driver for the memory-isolated ablation sweep.

Calls run_single_ablation.py via subprocess for each parameter value.
Each call is a fresh OS process -- if one crashes or OOMs, you get a
"CRASHED"/timeout row instead of losing the whole sweep (or your kernel).

Usage from a notebook:

    from ablation_driver import run_sweep

    results = run_sweep(
        npy_path="path/to/01_control_C2.npy",
        config_path="path/to/cell_config.json",
        vol_x=1600, vol_y=1200,
        n_frames=15,          # subset for faster iteration; set None for full data
        sweeps={
            "MotionModel.accuracy": [1.0, 2.0, 7.5, 15.0, 30.0, 50.0],
            "MotionModel.dt": [0.1, 0.5, 1.0, 2.0, 5.0],
            "MotionModel.max_lost": [2, 5, 10, 20],
        },
    )
    import pandas as pd
    pd.DataFrame(results)
"""

import subprocess
import json
import sys


WORKER_SCRIPT = "run_single_ablation.py"  # adjust path if not colocated


def run_one(config_path, vol_x, vol_y, objects_pkl_path=None, npy_path=None,
            df_path=None, well_id=None, param=None, value=None, n_frames=None,
            features=(), tracking_updates=("MOTION",), crop=None,
            update_method="EXACT", search_radius=None, optimize=False, timeout=300):
    cmd = [
        sys.executable, WORKER_SCRIPT,
        "--config_path", config_path,
        "--vol_x", str(vol_x),
        "--vol_y", str(vol_y),
        "--tracking_updates", ",".join(tracking_updates),
        "--update_method", update_method,
    ]
    if optimize:
        cmd += ["--optimize"]
    if search_radius is not None:
        cmd += ["--search_radius", str(search_radius)]
    if crop is not None:
        cmd += ["--crop", ",".join(str(v) for v in crop)]
    if objects_pkl_path is not None:
        cmd += ["--objects_pkl_path", objects_pkl_path]
    elif df_path is not None:
        cmd += ["--df_path", df_path]
        if well_id is not None:
            cmd += ["--well_id", str(well_id)]
    elif npy_path is not None:
        cmd += ["--npy_path", npy_path]
    else:
        raise ValueError("Must supply objects_pkl_path (fastest), df_path (fast), "
                          "or npy_path (slow)")

    if param is not None:
        cmd += ["--param", param, "--value", json.dumps(value)]
    if n_frames is not None:
        cmd += ["--n_frames", str(n_frames)]
    if features:
        cmd += ["--features", ",".join(features)]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT: {param}={value}")
        return {"param": param, "value": value, "status": "timeout"}

    if proc.returncode != 0:
        print(f"  CRASHED (exit {proc.returncode}): {param}={value}")
        # last bit of stderr often has the real error (e.g. segfault message,
        # OOM killer message won't show here but a nonzero/None returncode
        # with empty stdout usually means OOM-kill)
        print("   stderr tail:", proc.stderr[-500:] if proc.stderr else "(empty)")
        return {"param": param, "value": value, "status": "crashed",
                "returncode": proc.returncode}

    result_line = None
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT_JSON:"):
            result_line = line[len("RESULT_JSON:"):]
    if result_line is None:
        print(f"  NO RESULT parsed for {param}={value}")
        print("   stdout tail:", proc.stdout[-500:])
        return {"param": param, "value": value, "status": "no_result"}

    result = json.loads(result_line)
    result["status"] = "ok"
    debug_lines = [l for l in proc.stderr.splitlines() if l.startswith("DEBUG:")]
    for l in debug_lines:
        print(f"   {l}")
    print(f"  {param}={value} -> n_tracks={result['n_tracks']} "
          f"mean_len={result['mean_len']:.2f} median_len={result['median_len']:.1f} "
          f"frac_len1={result['frac_len1']:.2%} frac_len_lt5={result['frac_len_lt5']:.2%} "
          f"frac_len_ge20={result['frac_len_ge20']:.2%} max_len={result['max_len']}")
    return result


def objects_to_dicts(objects, feature_names=()):
    """
    Convert a list of PyTrackObject (NOT picklable -- ctypes pointers)
    into a list of plain dicts (picklable). Call this once on objects
    you've already built with segmentation_to_objects, pickle the result,
    and pass that path as objects_pkl_path to run_sweep/run_one.
    """
    dicts = []
    for o in objects:
        d = {"ID": int(o.ID), "t": int(o.t), "x": float(o.x), "y": float(o.y),
             "z": float(o.z), "label": int(o.label)}
        for feat in feature_names:
            d[feat] = float(getattr(o, feat))
        dicts.append(d)
    return dicts


def save_objects_pkl(objects, path, feature_names=()):
    """Convenience: objects_to_dicts + pickle.dump in one call."""
    import pickle
    dicts = objects_to_dicts(objects, feature_names)
    with open(path, "wb") as f:
        pickle.dump(dicts, f)
    return path


def run_sweep(config_path, vol_x, vol_y, sweeps: dict,
              objects_pkl_path=None, npy_path=None, df_path=None, well_id=None,
              n_frames=None, features=(), tracking_updates=("MOTION",),
              crop=None, update_method="EXACT", search_radius=None,
              optimize=False, timeout=300):
    """
    sweeps: dict of {dotted_param_name: [values...]}
    Also runs one baseline (unmodified config) call first.
    Returns a flat list of result dicts -- wrap in pd.DataFrame(...) after.
    """
    all_results = []

    print("Baseline (unmodified config):")
    all_results.append(
        run_one(config_path, vol_x, vol_y, objects_pkl_path=objects_pkl_path,
                npy_path=npy_path, df_path=df_path, well_id=well_id,
                param="baseline", value=None, n_frames=n_frames, features=features,
                tracking_updates=tracking_updates, crop=crop,
                update_method=update_method, search_radius=search_radius,
                optimize=optimize, timeout=timeout)
    )

    for param, values in sweeps.items():
        print(f"\n{param} sweep:")
        for v in values:
            all_results.append(
                run_one(config_path, vol_x, vol_y, objects_pkl_path=objects_pkl_path,
                        npy_path=npy_path, df_path=df_path, well_id=well_id,
                        param=param, value=v, n_frames=n_frames, features=features,
                        tracking_updates=tracking_updates, crop=crop,
                        update_method=update_method, search_radius=search_radius,
                        optimize=optimize, timeout=timeout)
            )

    return all_results