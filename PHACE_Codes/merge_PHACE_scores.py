import sys
import os
import numpy as np
import pandas as pd
import pyreadr

def _maybe_int32(arr: np.ndarray) -> np.ndarray:
    """
    Cast to int32 if values are integer-like; otherwise return original array.
    This helps shrink the .npz size for position columns while staying safe.
    """
    if arr.size == 0:
        return arr.astype(np.int32, copy=False)
    rounded = np.rint(arr)
    if np.allclose(arr, rounded, rtol=0.0, atol=0.0, equal_nan=True):
        # Guard against values outside int32 range
        if np.nanmin(rounded) >= np.iinfo(np.int32).min and np.nanmax(rounded) <= np.iinfo(np.int32).max:
            return rounded.astype(np.int32, copy=False)
    return arr

def main():
    if len(sys.argv) not in (4, 5):
        print(f"Usage: python {sys.argv[0]} <id> <num_parts> <output_base> [csv|--csv|<csv_output_path>]")
        sys.exit(1)

    id = sys.argv[1]
    num_parts = int(sys.argv[2])
    output_base = sys.argv[3]
    csv_arg = sys.argv[4] if len(sys.argv) == 5 else None

    if output_base.endswith(".feather"):
        output_base = output_base[:-8]
    elif output_base.endswith(".npz"):
        output_base = output_base[:-4]

    write_csv = False
    csv_file = f"{output_base}.csv"
    if csv_arg is not None:
        normalized = csv_arg.strip().lower()
        if normalized in {"csv", "--csv", "-c", "true", "1", "yes", "y"}:
            write_csv = True
        elif normalized in {"false", "0", "no", "n"}:
            write_csv = False
        else:
            write_csv = True
            csv_file = csv_arg if csv_arg.endswith(".csv") else f"{csv_arg}.csv"

    all_mats = []

    print(f"Merging {num_parts} parts for ID: {id}")
    
    for i in range(1, num_parts + 1):
        file_name = f"PHACE_scores/{id}/{id}_PHACE_part{i}.RData"
        if os.path.exists(file_name):
            try:
                # Read RData file
                result = pyreadr.read_r(file_name)
                # Get the first (or only) dataframe from the RData file
                # You may need to adjust the key if your RData has a specific object name
                mat = list(result.values())[0].values
                all_mats.append(mat)
                print(f"Loaded part {i}/{num_parts} (found {len(mat)} pairs)")
            except Exception as e:
                print(f"Could not load {file_name}: {e}")
        else:
            print(f"Warning: Part {i} not found ({file_name})")

    if not all_mats:
        print("No parts found. Exiting.")
        sys.exit(1)

    # Equivalent to R's do.call(rbind, ...)
    final_mat_np = np.vstack(all_mats)
    
    # Use pandas for the column-based correction logic
    final_mat = pd.DataFrame(final_mat_np, columns=["Pos1", "Pos2", "PHACE_Score"])
    
    # 'mean_score' calculation (same as before)
    total_pos = int(final_mat['Pos2'].max())
    ones_to_add = np.ones(total_pos)
    scores_with_ones = np.concatenate([final_mat['PHACE_Score'].values, ones_to_add])
    
    mean_score = np.mean(scores_with_ones)
    print(f"Calculated mean score (with {total_pos} ones): {mean_score}")

    # The 'PHACE_Corrected' column
    final_mat['PHACE_Corrected'] = final_mat['PHACE_Score'] - mean_score

    # Save the final result in Feather format
    feather_file = f"{output_base}.feather"
    final_mat.to_feather(feather_file)
    print(f"Saved Feather binary to {feather_file}")

    # Also save a compact NumPy .npz (compressed) for smaller storage
    npz_file = f"{output_base}.npz"
    pos1 = _maybe_int32(final_mat["Pos1"].to_numpy())
    pos2 = _maybe_int32(final_mat["Pos2"].to_numpy())
    score = final_mat["PHACE_Score"].to_numpy()
    corrected = final_mat["PHACE_Corrected"].to_numpy()
    np.savez_compressed(
        npz_file,
        Pos1=pos1,
        Pos2=pos2,
        PHACE_Score=score,
        PHACE_Corrected=corrected,
        mean_score=np.array(mean_score),
    )
    print(f"Saved compressed NumPy archive to {npz_file}")

    if write_csv:
        final_mat.to_csv(csv_file, index=False)
        print(f"Saved CSV to {csv_file}")

    print("Merging completed.")

if __name__ == "__main__":
    main()
