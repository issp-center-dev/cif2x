#!/usr/bin/env python
import os
import sys
import toml
import argparse
import pymatgen as mg
from pymatgen.ext.matproj import MPRester

def load_input_file(filename):
    """
    Load and validate the input TOML file.
    """
    if not os.path.exists(filename):
        print(f"Error: {filename} is not found.")
        sys.exit(1)

    info = toml.load(filename)
    for key in ["user_id", "criteria", "properties"]:
        if key not in info:
            print(f"Error: {key} is not found in {filename}.")
            sys.exit(1)
    return info

def save_properties(data, dir_name, properties):
    """
    Save the specified properties of a material into separate files.
    Exclude 'material_id' from the properties to be saved.
    """
    for prop in properties:
        with open(os.path.join(dir_name, prop), mode="w") as fp:
            fp.write(str(data[prop]) + "\n")

def main():
    parser = argparse.ArgumentParser(description='Get CIF files from Materials Project')
    parser.add_argument('-i', '--input', default="input.toml", help='input file name')
    args = parser.parse_args()

    info = load_input_file(args.input)

    user_id = info["user_id"]
    criteria = info["criteria"]
    properties = info["properties"]
    print(info)
    if "material_id" not in properties:
        properties.append("material_id")

    with MPRester(user_id) as m:
        data = m.query(criteria=criteria, properties=properties)
        count = 0
        for c in data:
            mpid = c["material_id"]
            if "mvc" in mpid:
                continue
            dir_name = os.path.join("mpall{}".format(mpid[3]), mpid)
            os.makedirs(dir_name, exist_ok=True)
            count += 1
            try:
                save_properties(c, dir_name, properties) # Save properties
                with open(os.path.join(dir_name, c["pretty_formula"] + ".cif"), mode="w") as fp:
                    fp.write(c["cif"])
            except Exception as e:
                print(f"Error with {c['material_id']}: {e}")

        print(count)

if __name__ == "__main__":
    main()
