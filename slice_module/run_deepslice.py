import os
import argparse
import traceback

# Attempt to import DeepSlice and provide a helpful error if it fails.
try:
    from DeepSlice import DSModel
except ImportError:
    print("ERROR: DeepSlice module not found.")
    print("Please ensure that:")
    print("1. The 'deepslice' conda environment is activated when you run this script OR that MATLAB is configured to use this environment.")
    print("2. DeepSlice is correctly installed in that environment.")
    print("3. If running this script directly, you might need to be in a directory where Python can find the DeepSlice module, or have the DeepSlice parent directory in your PYTHONPATH.")
    exit(1)

def process_folder(folder_path_abs, species_name, use_ensemble, use_section_numbers, deepslice_parent_dir_abs):
    """
    Processes a folder of images using DeepSlice.

    Args:
        folder_path_abs (str): Absolute path to the image folder.
        species_name (str): Species for the model ('mouse' or 'rat').
        use_ensemble (bool): Whether to use the ensemble model for prediction.
        use_section_numbers (bool): Whether section numbers are included in filenames.
        deepslice_parent_dir_abs (str): Absolute path to the DeepSlice parent directory.
    """
    original_dir = os.getcwd()
    print(f"--- Python Script Start ---")
    print(f"Original working directory: {original_dir}")
    print(f"Received image folder path: {folder_path_abs}")
    print(f"Received DeepSlice parent directory: {deepslice_parent_dir_abs}")
    print(f"Species: {species_name}, Ensemble: {use_ensemble}, Section Numbers: {use_section_numbers}")

    if not os.path.isdir(deepslice_parent_dir_abs):
        print(f"ERROR: DeepSlice parent directory not found: {deepslice_parent_dir_abs}")
        return

    try:
        # Change to the DeepSlice parent directory.
        # This is often crucial if DeepSlice relies on relative paths for its internal operations or model loading,
        # as suggested by the notebook's os.chdir('../../').
        os.chdir(deepslice_parent_dir_abs)
        print(f"Changed working directory to DeepSlice parent: {os.getcwd()}")

        # Initialize DeepSlice Model
        print(f"Initializing DSModel for species: {species_name}...")
        Model = DSModel(species_name)
        print("DSModel initialized.")

        # Predict
        # The folder_path_abs is already absolute, so it's fine to use directly.
        print(f"Starting prediction for folder: {folder_path_abs}...")
        Model.predict(folder_path_abs, ensemble=use_ensemble, section_numbers=use_section_numbers)
        print("Prediction complete.")

        # Propagate angles
        print("Propagating angles...")
        Model.propagate_angles()
        print("Angle propagation complete.")

        # Enforce index spacing
        # The notebook mentions options like Model.enforce_index_order() or
        # Model.enforce_index_spacing(section_thickness=25).
        # Here, we use the default enforce_index_spacing() as per the primary example in the notebook.
        print("Enforcing index spacing...")
        Model.enforce_index_spacing()
        print("Index spacing enforced.")

        # Save predictions
        # The notebook saves to folderpath + "results".
        # We'll create a slightly different name to distinguish these results.
        results_basename = os.path.join(folder_path_abs, "results_from_matlab_script")
        print(f"Saving predictions with base name: {results_basename}...")
        Model.save_predictions(results_basename) # DeepSlice adds .json and .xml extensions
        print(f"Predictions saved. Check for {results_basename}.json and {results_basename}.xml.")

        print("DeepSlice processing complete.")

    except Exception as e:
        print(f"ERROR: An error occurred during DeepSlice processing:")
        print(str(e))
        traceback.print_exc() # Prints the full traceback
    finally:
        # Crucial: change back to the original directory
        os.chdir(original_dir)
        print(f"Restored working directory to: {original_dir}")
        print(f"--- Python Script End ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DeepSlice processing. Called by MATLAB script.")
    parser.add_argument("folderpath", type=str, help="Absolute path to the image folder.")
    parser.add_argument("deepslice_dir", type=str, help="Absolute path to the DeepSlice parent directory (the one containing the 'DeepSlice' module/package).")
    parser.add_argument("--species", type=str, default="mouse", choices=["mouse", "rat"], help="Species for the model (default: mouse).")
    
    # For boolean flags, store_true makes them true if present, false otherwise.
    # We set defaults to True to match the notebook's direct usage.
    parser.add_argument("--ensemble", dest='ensemble', action='store_true', help="Use ensemble model (default: enabled).")
    parser.add_argument("--no-ensemble", dest='ensemble', action='store_false', help="Do NOT use ensemble model.")
    parser.set_defaults(ensemble=True)

    parser.add_argument("--section_numbers", dest='section_numbers', action='store_true', help="Indicates section numbers are in filenames (e.g., _sXXX) (default: enabled).")
    parser.add_argument("--no-section_numbers", dest='section_numbers', action='store_false', help="Indicates section numbers are NOT in filenames.")
    parser.set_defaults(section_numbers=True)

    args = parser.parse_args()

    # Ensure paths are absolute for robustness, though MATLAB should send absolute paths.
    abs_folder_path = os.path.abspath(args.folderpath)
    abs_deepslice_dir = os.path.abspath(args.deepslice_dir)

    if not os.path.isdir(abs_folder_path):
        print(f"CRITICAL ERROR: Image folder not found: {abs_folder_path}")
        exit(1)

    if not os.path.isdir(abs_deepslice_dir):
        print(f"CRITICAL ERROR: DeepSlice parent directory not found: {abs_deepslice_dir}")
        exit(1)
    
    # Call the main processing function
    process_folder(abs_folder_path, args.species, args.ensemble, args.section_numbers, abs_deepslice_dir)
