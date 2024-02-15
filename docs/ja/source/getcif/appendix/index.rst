.. _Ch:Appendix:

.. raw:: latex

   \appendix

================================================================
パラメータリスト
================================================================

検索条件 (properties)
----------------------------------------------------------------

properties に指定できる項目と、その項目がどのような値を取るかを以下にまとめます。

Materials Project API のクライアントアプリケーションの一つとして mp-api パッケージが Materials Project から公開されており、getcif はこのライブラリを利用してデータベースへの接続を行います。
以下は MPRester クラスの materials.summary.search メソッドのパラメータに対応します。
(以下の表は materials.summary.search のソースコードのコメントから転記しました。)

値の型の表記は次のとおりです。

- ``str``: 文字列型
- ``List[str]``:  文字列型のリスト
- ``str | List[str]``: 単一の文字列、または、文字列型のリスト
- ``int``: 整数型
- ``bool``: 真偽値 (true または false)
- ``Tuple[float,float]``: 実数値 2つからなる組 (リスト)
- ``Tuple[int,int]``: 整数値 2つからなる組 (リスト)
- ``CrystalSystem``: 結晶のタイプを表す文字列。Triclinic, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic のいずれか。
- ``List[HasProps]``: 特性値のタイプを表す文字列のリスト。特性値は emmet.core.summary に定義されている。以下のいずれかの値を取る。

    absorption,
    bandstructure,
    charge_density,
    chemenv,
    dielectric,
    dos,
    elasticity,
    electronic_structure,
    eos,
    grain_boundaries,
    insertion_electrodes,
    magnetism,
    materials,
    oxi_states,
    phonon,
    piezoelectric,
    provenance,
    substrates,
    surface_properties,
    thermo,
    xas

- ``Ordering``: 磁気秩序を表す文字列。FM, AFM, FiM, NM のいずれか。

値のリストは、YAML形式の箇条書きおよび ``[ ... ]`` にカンマ区切りで記述するほか、空白区切りで列挙する記法も可能です。

``Tuple`` で表される型は値の範囲 (min, max) の指定に使われます。値のリストとして記述するほか、空白区切りで ``min max`` のように記述することもできます。また、以下の表記も可能です。

     ``<= max``     : max 以下

     ``< max``      : max より小さい

     ``>= min``     : min 以上

     ``> min``      : min より大きい

     ``min ~ max``  : min 以上 max 以下

.. _getcif-cond-table:

.. list-table:: 検索条件のキーワード
    :widths: 30 20 60
    :header-rows: 1

    * - Keyword
      - Type
      - Description
    * - band_gap
      - Tuple[float,float]
      - Minimum and maximum band gap in eV to consider.
    * - chemsys
      - str | List[str]
      - A chemical system, list of chemical systems (e.g., Li-Fe-O, Si-\*, [Si-O, Li-Fe-P]), or single formula (e.g., Fe2O3, Si\*).
    * - crystal_system
      - CrystalSystem
      - Crystal system of material.
    * - density
      - Tuple[float,float]
      - Minimum and maximum density to consider.
    * - deprecated
      - bool
      - Whether the material is tagged as deprecated.
    * - e_electronic
      - Tuple[float,float]
      - Minimum and maximum electronic dielectric constant to consider.
    * - e_ionic
      - Tuple[float,float]
      - Minimum and maximum ionic dielectric constant to consider.
    * - e_total
      - Tuple[float,float]
      - Minimum and maximum total dielectric constant to consider.
    * - efermi
      - Tuple[float,float]
      - Minimum and maximum fermi energy in eV to consider.
    * - elastic_anisotropy
      - Tuple[float,float]
      - Minimum and maximum value to consider for the elastic anisotropy.
    * - elements
      - List[str]
      - A list of elements.
    * - energy_above_hull
      - Tuple[int,int]
      - Minimum and maximum energy above the hull in eV/atom to consider.
    * - equilibrium_reaction_energy
      - Tuple[float,float]
      - Minimum and maximum equilibrium reaction energy in eV/atom to consider.
    * - exclude_elements
      - List[str]
      - List of elements to exclude.
    * - formation_energy
      - Tuple[int,int]
      - Minimum and maximum formation energy in eV/atom to consider.
    * - formula
      - str | List[str]
      - A formula including anonymized formula or wild cards (e.g., Fe2O3, ABO3, Si\*). A list of chemical formulas can also be passed (e.g., [Fe2O3, ABO3]).
    * - g_reuss
      - Tuple[float,float]
      - Minimum and maximum value in GPa to consider for the Reuss average of the shear modulus.
    * - g_voigt
      - Tuple[float,float]
      - Minimum and maximum value in GPa to consider for the Voigt average of the shear modulus.
    * - g_vrh
      - Tuple[float,float]
      - Minimum and maximum value in GPa to consider for the Voigt-Reuss-Hill average of the shear modulus.
    * - has_props
      - List[HasProps]
      - The calculated properties available for the material.
    * - has_reconstructed
      - bool
      - Whether the entry has any reconstructed surfaces.
    * - is_gap_direct
      - bool
      - Whether the material has a direct band gap.
    * - is_metal
      - bool
      - Whether the material is considered a metal.
    * - is_stable
      - bool
      - Whether the material lies on the convex energy hull.
    * - k_reuss
      - Tuple[float,float]
      - Minimum and maximum value in GPa to consider for the Reuss average of the bulk modulus.
    * - k_voigt
      - Tuple[float,float]
      - Minimum and maximum value in GPa to consider for the Voigt average of the bulk modulus.
    * - k_vrh
      - Tuple[float,float]
      - Minimum and maximum value in GPa to consider for the Voigt-Reuss-Hill average of the bulk modulus.
    * - magnetic_ordering
      - Ordering
      - Magnetic ordering of the material.
    * - material_ids
      - List[str]
      - List of Materials Project IDs to return data for.
    * - n
      - Tuple[float,float]
      - Minimum and maximum refractive index to consider.
    * - num_elements
      - Tuple[int,int]
      - Minimum and maximum number of elements to consider.
    * - num_sites
      - Tuple[int,int]
      - Minimum and maximum number of sites to consider.
    * - num_magnetic_sites
      - Tuple[int,int]
      - Minimum and maximum number of magnetic sites to consider.
    * - num_unique_magnetic_sites
      - Tuple[int,int]
      - Minimum and maximum number of unique magnetic sites to consider.
    * - piezoelectric_modulus
      - Tuple[float,float]
      - Minimum and maximum piezoelectric modulus to consider.
    * - poisson_ratio
      - Tuple[float,float]
      - Minimum and maximum value to consider for Poisson's ratio.
    * - possible_species
      - List[str]
      - List of element symbols appended with oxidation states. (e.g. Cr2+,O2-)
    * - shape_factor
      - Tuple[float,float]
      - Minimum and maximum shape factor values to consider.
    * - spacegroup_number
      - int
      - Space group number of material.
    * - spacegroup_symbol
      - str
      - Space group symbol of the material in international short symbol notation.
    * - surface_energy_anisotropy
      - Tuple[float,float]
      - Minimum and maximum surface energy anisotropy values to consider.
    * - theoretical
      - bool
      - Whether the material is theoretical.
    * - total_energy
      - Tuple[int,int]
      - Minimum and maximum corrected total energy in eV/atom to consider.
    * - total_magnetization
      - Tuple[float,float]
      - Minimum and maximum total magnetization values to consider.
    * - total_magnetization_normalized_formula_units
      - Tuple[float,float]
      - Minimum and maximum total magnetization values normalized by formula units to consider.
    * - total_magnetization_normalized_vol
      - Tuple[float,float]
      - Minimum and maximum total magnetization values normalized by volume to consider.
    * - uncorrected_energy
      - Tuple[int,int]
      - Minimum and maximum uncorrected total energy in eV/atom to consider.
    * - volume
      - Tuple[float,float]
      - Minimum and maximum volume to consider.
    * - weighted_surface_energy
      - Tuple[float,float]
      - Minimum and maximum weighted surface energy in J/:math:`m^2` to consider.
    * - weighted_work_function
      - Tuple[float,float]
      - Minimum and maximum weighted work function in eV to consider.

..
.. .. list-table:: Unsupported search criteria for the properties section
..     :widths: 30 20 60
..    :header-rows: 1
..
..    * - Keyword
..      - Type
..      - Description
..    * - num_chunks
..      - int
..      - Maximum number of chunks of data to yield. None will yield all possible.
..    * - chunk_size
..      - int
..      - Number of data entries per chunk.
..    * - all_fields
..      - bool
..      - Whether to return all fields in the document. Defaults to True.
..    * - fields
..      - List[str]
..      - List of fields in SearchDoc to return data for. Default is material_id if all_fields is False.
..


出力項目 (fields)
----------------------------------------------------------------

fields に指定できる項目を以下に列挙します。

.. code:: pre

    band_gap
    bandstructure
    builder_meta
    bulk_modulus
    cbm
    chemsys
    composition
    composition_reduced
    database_IDs
    decomposes_to
    density
    density_atomic
    deprecated
    deprecation_reasons
    dos
    dos_energy_down
    dos_energy_up
    e_electronic
    e_ij_max
    e_ionic
    e_total
    efermi
    elements
    energy_above_hull
    energy_per_atom
    equilibrium_reaction_energy_per_atom
    es_source_calc_id
    formation_energy_per_atom
    formula_anonymous
    formula_pretty
    grain_boundaries
    has_props
    has_reconstructed
    homogeneous_poisson
    is_gap_direct
    is_magnetic
    is_metal
    is_stable
    last_updated
    material_id
    n
    nelements
    nsites
    num_magnetic_sites
    num_unique_magnetic_sites
    ordering
    origins
    possible_species
    property_name
    shape_factor
    shear_modulus
    structure
    surface_anisotropy
    symmetry
    task_ids
    theoretical
    total_magnetization
    total_magnetization_normalized_formula_units
    total_magnetization_normalized_vol
    types_of_magnetic_species
    uncorrected_energy_per_atom
    universal_anisotropy
    vbm
    volume
    warnings
    weighted_surface_energy
    weighted_surface_energy_EV_PER_ANG2
    weighted_work_function
    xas
