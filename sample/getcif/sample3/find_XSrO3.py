import os,sys
from ruamel.yaml import YAML
from getcif import QueryMaterialsProject
import copy

import logging
logger = logging.getLogger("main")

if len(sys.argv) < 2:
    print("usage: {} input_file".format(sys.argv[0]))
    sys.exit(0)
else:
    input_file = sys.argv[1]

yaml = YAML(typ="safe")
with open(input_file, "r") as fp:
    info = yaml.load(fp)

data = []
for elem in ["Co", "Ti", "V"]:
    info["properties"].update({"elements": [ elem, "Sr" ]})
    data += QueryMaterialsProject(info).run()

if len(data) > 0:
    print([[idx+1, v["material_id"], v["formula"]] for idx, v in enumerate(data)])
    
print("done.")
