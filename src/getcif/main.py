#!/usr/bin/env python3

import os
import sys
import copy
import io
import csv
import json
import argparse
import logging

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from ruamel.yaml import YAML, YAMLError
from mp_api.client import MPRester

# enum parameters
from emmet.core.symmetry import CrystalSystem
from emmet.core.summary import HasProps
from pymatgen.analysis.magnetism import Ordering

logger = logging.getLogger("getcif")

try:
    __version__ = version("htp-tools-getcif")
except PackageNotFoundError:
    __version__ = "develop"

class QueryMaterialsProject:
    query_table = {
        "band_gap": "tuple[float,float]",
        "chemsys": "str|list[str]",
        "crystal_system": "CrystalSystem",
        "density": "tuple[float,float]",
        "deprecated": "bool",
        "e_electronic": "tuple[float,float]",
        "e_ionic": "tuple[float,float]",
        "e_total": "tuple[float,float]",
        "efermi": "tuple[float,float]",
        "elastic_anisotropy": "tuple[float,float]",
        "elements": "list[str]",
        "energy_above_hull": "tuple[float,float]",
        "equilibrium_reaction_energy": "tuple[float,float]",
        "exclude_elements": "list[str]",
        "formation_energy": "tuple[float,float]",
        "formula": "str|list[str]",
        "g_reuss": "tuple[float,float]",
        "g_voigt": "tuple[float,float]",
        "g_vrh": "tuple[float,float]",
        "has_props": "list[HasProps]",
        "has_reconstructed": "bool",
        "is_gap_direct": "bool",
        "is_metal": "bool",
        "is_stable": "bool",
        "k_reuss": "tuple[float,float]",
        "k_voigt": "tuple[float,float]",
        "k_vrh": "tuple[float,float]",
        "magnetic_ordering": "Ordering",
        "material_ids": "list[str]",
        "n": "tuple[float,float]",
        "num_elements": "tuple[int,int]",
        "num_sites": "tuple[int,int]",
        "num_magnetic_sites": "tuple[int,int]",
        "num_unique_magnetic_sites": "tuple[int,int]",
        "piezoelectric_modulus": "tuple[float,float]",
        "poisson_ratio": "tuple[float,float]",
        "possible_species": "list[str]",
        "shape_factor": "tuple[float,float]",
        "spacegroup_number": "int",
        "spacegroup_symbol": "str",
        "surface_energy_anisotropy": "tuple[float,float]",
        "theoretical": "bool",
        "total_energy": "tuple[float,float]",
        "total_magnetization": "tuple[float,float]",
        "total_magnetization_normalized_formula_units": "tuple[float,float]",
        "total_magnetization_normalized_vol": "tuple[float,float]",
        "uncorrected_energy": "tuple[float,float]",
        "volume": "tuple[float,float]",
        "weighted_surface_energy": "tuple[float,float]",
        "weighted_work_function": "tuple[float,float]",
    }

    # mpr.materials.summary.available_fields
    available_fields = [
        "band_gap",
        "bandstructure",
        "builder_meta",
        "bulk_modulus",
        "cbm",
        "chemsys",
        "composition",
        "composition_reduced",
        "database_IDs",
        "decomposes_to",
        "density",
        "density_atomic",
        "deprecated",
        "deprecation_reasons",
        "dos",
        "dos_energy_down",
        "dos_energy_up",
        "e_electronic",
        "e_ij_max",
        "e_ionic",
        "e_total",
        "efermi",
        "elements",
        "energy_above_hull",
        "energy_per_atom",
        "equilibrium_reaction_energy_per_atom",
        "es_source_calc_id",
        "formation_energy_per_atom",
        "formula_anonymous",
        "formula_pretty",
        "grain_boundaries",
        "has_props",
        "has_reconstructed",
        "homogeneous_poisson",
        "is_gap_direct",
        "is_magnetic",
        "is_metal",
        "is_stable",
        "last_updated",
        "material_id",
        "n",
        "nelements",
        "nsites",
        "num_magnetic_sites",
        "num_unique_magnetic_sites",
        "ordering",
        "origins",
        "possible_species",
        "property_name",
        "shape_factor",
        "shear_modulus",
        "structure",
        "surface_anisotropy",
        "symmetry",
        "task_ids",
        "theoretical",
        "total_magnetization",
        "total_magnetization_normalized_formula_units",
        "total_magnetization_normalized_vol",
        "types_of_magnetic_species",
        "uncorrected_energy_per_atom",
        "universal_anisotropy",
        "vbm",
        "volume",
        "warnings",
        "weighted_surface_energy",
        "weighted_surface_energy_EV_PER_ANG2",
        "weighted_work_function",
        "xas",
    ]

    def __init__(self, info):
        self.info = copy.deepcopy(info)
        self._setup_dbinfo(self.info.get("database", {}))
        self._setup_option(self.info.get("option", {}))

    def _setup_dbinfo(self, info):
        self.dbinfo = info

        # setup api key
        # 1. read from api_key_file (default "materials_project.key") if exists
        # 2. taken from environment variable or pymetgen settings (leave api_key None)

        api_key = None

        api_key_file = info.get("api_key_file", "materials_project.key")
        if api_key_file.endswith(".key") and Path(api_key_file).exists():
            with open(Path(api_key_file), "r", encoding="utf-8") as fp:
                data = [s.strip() for s in fp.readlines() if not s.strip().startswith("#")]
                if data:
                    api_key = data[0]
        if not api_key:
            logger.debug("api_key not set. use environment variable or pymatgen settings")

        self.api_key = api_key

    def _setup_option(self, info):
        self.output_dir = info.get("output_dir", "")
        self.dry_run = info.get("dry_run", False)
        # symprec: default value 0.1 used in Materials Project to determine symmetry
        self.symprec = info.get("symprec", 0.1)
        if self.symprec == 0:
            self.symprec = None

    def _find_query(self, info):
        props = self._find_properties(info.get("properties", {}))

        fields = self._find_fields(info.get("fields", ""))

        props.update({"fields": fields})
        return props

    def _find_fields(self, info):
        if isinstance(info, str):
            fields = info.split()
        elif isinstance(info, list):
            fields = info
        else:
            raise ValueError("invalid fields parameter")

        if not "material_id" in fields:
            fields.append("material_id")
        if not "formula_pretty" in fields:
            fields.append("formula_pretty")

        # check
        err = 0
        for field in fields:
            if not field in QueryMaterialsProject.available_fields:
                logger.error("unknown field name {}".format(field))
                err += 1

        if err > 0:
            logger.error("fields check failed")
            raise ValueError("fields check failed")

        logger.info("query: fields={}".format(fields))
        return fields

    def _find_properties(self, info):
        props = info

        def _find_val_or_none(val, typ=float):
            if isinstance(val, str):
                if val.lower() == "none":
                    return None
            return typ(val)

        def _find_range(val, typ=float):
            """
            accepted patterns
            keyword: < 1.0
            keyword: <= 1.0
            keyword: > 0.5
            keyword: >= 0.5
            keyword: 0.5 ~ 1.0
            keyword: 0.5 1.0
            keyword: = 0.5
            not accepted
            keyword: 1.0 <  #NG
            """
            if any(symbol in val for symbol in ("<",">","=","~")):
                w = [s.strip() for s in val.split()]
                if len(w) == 2:
                    if w[0] == "<=":
                        return (None, typ(w[1]))
                    if w[0] == "<":
                        if typ == int:
                            return (None, typ(w[1])-1)
                        else:
                            return (None, typ(w[1]))
                    if w[0] == ">=":
                        return (typ(w[1]), None)
                    if w[0] == ">":
                        if typ == int:
                            return (typ(w[1])+1, None)
                        else:
                            return (typ(w[1]), None)
                    if w[0] == "=":
                        return (typ(w[1]), typ(w[1]))
                elif len(w) == 3:
                    if w[1] == "~":
                        return (typ(w[0]), typ(w[2]))
            else:
                w = [typ(s) for s in val.split()]
                if len(w) == 2:
                    return tuple(w[0:2])
            raise ValueError("illegal string: {}".format(val))

        # format and check
        err = 0
        for prop,value in props.items():
            typ = QueryMaterialsProject.query_table.get(prop, None)
            if typ is None:
                logger.error("unknown query key {}".format(prop))
                err += 1
            elif typ == "bool":
                pass
            elif typ == "int":
                pass
            elif typ == "str":
                pass
            elif typ == "list[str]":
                if isinstance(value, str):
                    props[prop] = value.split()
            elif typ == "str|list[str]":
                if value == "":
                    pass
                else:
                    v = value.split()
                    props[prop] = v if len(v) > 1 else v[0]
            elif typ == "tuple[int,int]":
                if isinstance(value, list):
                    props[prop] = tuple(_find_val_or_none(s, int) for s in value[0:2])
                elif isinstance(value, str):
                    props[prop] = _find_range(value, int)
            elif typ == "tuple[float,float]":
                if isinstance(value, list):
                    props[prop] = tuple(_find_val_or_none(s) for s in value[0:2])
                elif isinstance(value, str):
                    props[prop] = _find_range(value)
            elif typ == "list[HasProps]":  # for has_props
                if isinstance(value, list):
                    props[prop] = [HasProps[v.lower()] for v in value]
                elif isinstance(value, str):
                    props[prop] = [HasProps[v.lower()] for v in value.split()]
            elif typ == "CrystalSystem":  # for crystal_system
                props[prop] = CrystalSystem[value.capitalize()]
            elif typ == "Ordering":  # for magnetic_ordering
                props[prop] = Ordering[value]
            else:
                logger.error("unknown query type {} for key {}".format(typ, prop))
            logger.debug("prop={}, value={}, type={} -> {}".format(prop,value,typ,props[prop]))

        if err > 0:
            logger.error("properties check failed")
            raise ValueError("properties check failed")

        logger.info("query: properties={}".format(props))
        return props

    def _do_query(self, query):
        if not self.dry_run:
            try:
                with MPRester(api_key=self.api_key, mute_progress_bars=True) as mpr:
                    docs = mpr.materials.summary.search(**query)
                    material_ids = [ doc.material_id for doc in docs ]
                    logger.info("result: number of entries={}".format(len(docs)))
                    logger.info("result: material_ids={}".format(material_ids))
            except Exception as e:
                logger.error("{}".format(e))
                raise
            return docs
        else:
            logger.info("dry run.")
            print(query)
            return []

    def _do_summary(self, docs, fields):
        results = []
        symprec = self.symprec

        for idx, doc in enumerate(docs):
            m_id = str(doc.material_id)
            m_formula = doc.formula_pretty
            d = dict(doc)

            data = { "material_id": m_id, "formula": m_formula }

            data_dir = Path(self.output_dir, m_id)
            os.makedirs(data_dir, exist_ok=True)

            for field in fields:
                if field == "material_id":
                    pass
                elif field == "structure":
                    if symprec is not None:
                        doc.structure.to(Path(data_dir, "structure.cif"), fmt="cif", symprec=symprec)
                    else:
                        doc.structure.to(Path(data_dir, "structure.cif"), fmt="cif")
                elif field == "formula_pretty":
                    with open(Path(data_dir, "formula"), "w", encoding="utf-8") as fp:
                        fp.write(str(d[field]) + "\n")
                else:
                    with open(Path(data_dir, field), "w", encoding="utf-8") as fp:
                        fp.write(str(d[field]) + "\n")
                    data.update({field: d[field]})

            # # export summary in json format
            # with open(Path(data_dir, "summary.json"), "w") as fp:
            #     fp.write(json.dumps(data))

            results.append(data)
        return results

    def _do_report(self, results, output_file=None, fmt="text"):
        if len(results) < 1:
            return None

        if fmt == "text":
            fields = results[0].keys()
            retv = "  ".join(fields) + "\n"
            for result in results:
                retv += "  ".join([str(result[field]) for field in fields]) + "\n"
        elif fmt == "csv":
            output = io.StringIO()

            writer = csv.DictWriter(output, fieldnames=results[0].keys())
            writer.writeheader()
            for result in results:
                writer.writerow(result)

            retv = output.getvalue()
        elif fmt == "json":
            retv = json.dumps(results)
        else:
            logger.error("unknown format {}".format(fmt))
            return None

        if output_file is None:
            print(retv)
        else:
            with open(output_file, "w", encoding="utf-8") as fp:
                fp.write(retv)
        return None

    def run(self, output_file=None, fmt="text"):
        query = self._find_query(self.info)
        docs = self._do_query(query)
        results = self._do_summary(docs, query["fields"])
        self._do_report(results, output_file, fmt)
        return results

    @classmethod
    def show_info(cls, file=None):
        if file is None:
            file = sys.stdout
        file.write("Query parameters:\n")
        file.write("  Properties:\n")
        for k, v in QueryMaterialsProject.query_table.items():
            file.write("    {:24s} {}\n".format(k, v))
        file.write("\n")
        file.write("  Fields:\n")
        for k in QueryMaterialsProject.available_fields:
            file.write("    {}\n".format(k))

def main():
    class ArgumentParser(argparse.ArgumentParser):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
        def print_help(self, file=None):
            if file is None:
                file = sys.stdout
            super().print_help(file)
            file.write("\n")
            QueryMaterialsProject.show_info(file)

    parser = ArgumentParser(prog="getcif")
    parser.add_argument("input_file", action="store",
                        help="input parameter file (input.yaml)")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="increase output verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0,
                        help="decrease output verbosity")
    parser.add_argument("--dry-run", action="store_true", default=False,
                        help="dry run")
    parser.add_argument("--version", action="version",
                        version="%(prog)s version {}".format(__version__))

    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING-(args.verbose-args.quiet)*10)

    try:
        yaml = YAML(typ="safe")
        with open(Path(args.input_file), mode="r", encoding="utf-8") as fp:
            info_dict = yaml.load(fp)
    except FileNotFoundError as e:
        logger.error("{}".format(e))
        sys.exit(1)
    except YAMLError as e:
        logger.error("{}".format(e))
        sys.exit(1)

    if info_dict is None:
        logger.error("input file is empty")
        # raise ValueError("empty input file")
        sys.exit(1)

    if args.dry_run:
        if "option" in info_dict:
            info_dict["option"].update({"dry_run": True})
        else:
            info_dict["option"] = {"dry_run": True}

    #try:
    stat = QueryMaterialsProject(info_dict).run()
    #except Exception as e:
    #    logger.error("{}".format(e))
    #    sys.exit(1)

if __name__ == "__main__":
    main()
