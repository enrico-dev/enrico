"""Script to be used to combine separate assembly VTK visualization
files into one single VTK file.
"""

import argparse
import warnings
import glob
from vtk import vtkAppendFilter, vtkUnstructuredGridWriter, \
    vtkUnstructuredGrid, vtkUnstructuredGridReader


def parse_args():
    parser = argparse.ArgumentParser(description="Combine separate assembly VTK " +\
        "visualization files into one single VTK file.")
    parser.add_argument('-f', '--files', type=str, nargs='+',
                        help="vtk files to combine")
    parser.add_argument('-o', '--output', type=str, help="name for output file")
    return parser.parse_args()


def check_input(files):
    """Check files are all type *.vtk"""
    ferror_list = []
    for f in files:
        if f[-4:] != ".vtk":
            ferror_list.append(f)
    if len(ferror_list) > 0:
        raise Exception("Files are not '*.vtk' type: {}".format(ferror_list))


def check_output(output):
    """Check that output name is valid type. Use default name if not."""
    if output[-4:] != ".vtk":
        warnings.warn("Output filename is not valid. " +\
            "Using default name: combined_output.vtk")
        output = 'combined_output.vtk'
    return output


def combine_files(files, outname):
    """combine all files into single vtk file named outname."""
    reader = vtkUnstructuredGridReader()
    append = vtkAppendFilter()
    for file in files:
        reader.SetFileName(file)
        reader.Update()
        unstructured = vtkUnstructuredGrid()
        unstructured.ShallowCopy(reader.GetOutput())
        append.AddInputData(unstructured)
    append.Update()
    writer = vtkUnstructuredGridWriter()
    writer.SetFileName(outname)
    writer.SetInputData(append.GetOutput())
    writer.Write()


if __name__ == "__main__":
    args = parse_args()

    # get file names
    if args.files:
        filenames = args.files
    else:
        filenames = glob.glob("*.vtk")

    # get output name
    if args.output:
        output = check_output(args.output)
    else:
        output = "combined_output.vtk"

    check_input(filenames)

    print("Combining files: {}".format(filenames))
    combine_files(filenames, output)