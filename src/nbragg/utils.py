import NCrystal as NC
from pathlib import Path

# Initialize materials as a global dictionary
materials = {}

def make_materials_dict():
    # populate the materials dictionary based on the available ncmat files in your system
    mat_dict = {}
    ncmat_files = NC.browseFiles()
    for mat in ncmat_files:
        if mat.name.endswith(".ncmat"):
            fullname = mat.name
            name = ""
            formula = ""
            space_group = ""
            splitname = fullname.replace(".ncmat", "").split("_")
            if len(splitname) == 3:
                formula, space_group, name = splitname
            elif len(splitname) == 1:
                formula = splitname[0]
                name = formula
            elif len(splitname) == 2:
                formula, space_group = splitname
                name = formula
            elif len(splitname) > 3:
                formula, space_group, *names = splitname
                name = "_".join(names)
            else:
                continue
            mat_dict[fullname] = {
                "name": name,
                "filename": fullname,
                "formula": formula,
                "space_group": space_group
            }
    return mat_dict

def initialize_materials():
    global materials
    materials.update(make_materials_dict())

def register_material(filename):
    filename = Path(filename)
    with open(filename, "r") as fid:
        NC.registerInMemoryFileData(filename.name, fid.read())
    initialize_materials()

# Initialize the materials dictionary at module import
initialize_materials()
