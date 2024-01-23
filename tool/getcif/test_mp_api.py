from pymatgen.ext.matproj import MPRester

# APIキーを設定（適切なキーに置き換えてください）
api_key = ""

# MPRester インスタンスの作成
with MPRester(api_key) as mpr:
    search_results = mpr.materials.search(material_ids=["mp-149", "mp-4163"], fields=["structure"])
    for i, doc in enumerate(search_results):
        structure = doc.structure
        cif_filename = f"material_{i}.cif"
        structure.to(fmt="cif", filename=cif_filename)
