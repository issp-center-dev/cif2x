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

値の型の表記は次のとおりです。

- ``str``: 文字列型
- ``list[str]``:  文字列型のリスト
- ``str|list[str]``: 単一の文字列、または、文字列型のリスト
- ``int``: 整数型
- ``bool``: 真偽値 (true または false)
- ``tuple[float,float]``: 実数値 2つからなる組 (リスト)
- ``tuple[int,int]``: 整数値 2つからなる組 (リスト)
- ``CrystalSystem``: 結晶のタイプを表す文字列。Triclinic, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic のいずれか。
- ``list[HasProps]``: 特性値のタイプを表す文字列のリスト。特性値は emmet.core.summary に定義されている。以下のいずれかの値を取る。

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

``tuple`` で表される型は値の範囲 (min, max) の指定に使われます。値のリストとして記述するほか、空白区切りで ``min max`` のように記述することもできます。また、以下の表記も可能です。

     ``< max``      : max 以下

     ``> min``      : min 以上

     ``min ~ max``  : min 以上 max 以下

.. code:: pre

    band_gap                 tuple[float,float]
    chemsys                  str|list[str]
    crystal_system           CrystalSystem
    density                  tuple[float,float]
    deprecated               bool
    e_electronic             tuple[float,float]
    e_ionic                  tuple[float,float]
    e_total                  tuple[float,float]
    efermi                   tuple[float,float]
    elastic_anisotropy       tuple[float,float]
    elements                 list[str]
    energy_above_hull        tuple[float,float]
    equilibrium_reaction_energy tuple[float,float]
    exclude_elements         list[str]
    formation_energy         tuple[float,float]
    formula                  str|list[str]
    g_reuss                  tuple[float,float]
    g_voigt                  tuple[float,float]
    g_vrh                    tuple[float,float]
    has_props                list[HasProps]
    has_reconstructed        bool
    is_gap_direct            bool
    is_metal                 bool
    is_stable                bool
    k_reuss                  tuple[float,float]
    k_voigt                  tuple[float,float]
    k_vrh                    tuple[float,float]
    magnetic_ordering        Ordering
    material_ids             list[str]
    n                        tuple[float,float]
    num_elements             tuple[int,int]
    num_sites                tuple[int,int]
    num_magnetic_sites       tuple[int,int]
    num_unique_magnetic_sites tuple[int,int]
    piezoelectric_modulus    tuple[float,float]
    poisson_ratio            tuple[float,float]
    possible_species         list[str]
    shape_factor             tuple[float,float]
    spacegroup_number        int
    spacegroup_symbol        str
    surface_energy_anisotropy tuple[float,float]
    theoretical              bool
    total_energy             tuple[float,float]
    total_magnetization      tuple[float,float]
    total_magnetization_normalized_formula_units tuple[float,float]
    total_magnetization_normalized_vol tuple[float,float]
    uncorrected_energy       tuple[float,float]
    volume                   tuple[float,float]
    weighted_surface_energy  tuple[float,float]
    weighted_work_function   tuple[float,float]


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
