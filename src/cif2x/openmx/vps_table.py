# VPS,  Valence electrons,  Quick,  Standard,  Precise
VPS_TABLE = [
  ["E",  "E",          0.0, "Kr10.0-s1p1",     "Kr10.0-s2p1d1",   "Kr10.0-s2p2d1f1",],
  ["H",  "H_PBE19",    1.0, "H5.0-s2",         "H6.0-s2p1",       "H7.0-s2p2d1",    ],
  ["He", "He_PBE19",   2.0, "He8.0-s1p1",      "He8.0-s2p1",      "He10.0-s2p2d1",  ],
  ["Li", "Li_PBE19",   3.0, "Li8.0-s3p1",      "Li8.0-s3p2",      "Li8.0-s3p2d1",   ],
  ["Be", "Be_PBE19",   2.0, "Be7.0-s2p1",      "Be7.0-s2p2",      "Be7.0-s3p2d1",   ],
  ["B",  "B_PBE19",    3.0, "B7.0-s2p2",       "B7.0-s2p2d1",     "B7.0-s3p2d2",    ],
  ["C",  "C_PBE19",    4.0, "C6.0-s2p2",       "C6.0-s2p2d1",     "C6.0-s3p2d2",    ],
  ["N",  "N_PBE19",    5.0, "N6.0-s2p2",       "N6.0-s2p2d1",     "N6.0-s3p2d2",    ],
  ["O",  "O_PBE19",    6.0, "O6.0-s2p2",       "O6.0-s2p2d1",     "O6.0-s3p2d2",    ],
  ["F",  "F_PBE19",    7.0, "F6.0-s2p2",       "F6.0-s2p2d1",     "F6.0-s3p3d2f1",  ],
  ["Ne", "Ne_PBE19",   8.0, "Ne9.0-s2p2",      "Ne9.0-s2p2d1",    "Ne9.0-s3p2d2",   ],
  ["Na", "Na_PBE19",   9.0, "Na9.0-s3p2",      "Na9.0-s3p2d1",    "Na9.0-s3p2d2",   ],
  ["Mg", "Mg_PBE19",   8.0, "Mg9.0-s2p2",      "Mg9.0-s3p2d1",    "Mg9.0-s3p2d2",   ],
  ["Al", "Al_PBE19",   3.0, "Al7.0-s2p1d1",    "Al7.0-s2p2d1",    "Al7.0-s3p2d2",   ],
  ["Si", "Si_PBE19",   4.0, "Si7.0-s2p1d1",    "Si7.0-s2p2d1",    "Si7.0-s3p3d2",   ],
  ["P",  "P_PBE19",    5.0, "P7.0-s2p2d1",     "P7.0-s2p2d1f1",   "P7.0-s3p2d2f1",  ],
  ["S",  "S_PBE19",    6.0, "S7.0-s2p2d1",     "S7.0-s2p2d1f1",   "S7.0-s3p2d2f1",  ],
  ["Cl", "Cl_PBE19",   7.0, "Cl7.0-s2p2d1",    "Cl7.0-s2p2d1f1",  "Cl7.0-s3p2d2f1", ],
  ["Ar", "Ar_PBE19",   8.0, "Ar9.0-s2p2d1",    "Ar9.0-s2p2d1f1",  "Ar9.0-s3p2d2f1", ],
  ["K",  "K_PBE19",    9.0, "K10.0-s3p2",      "K10.0-s3p2d1",    "K10.0-s3p2d2",   ],
  ["Ca", "Ca_PBE19",  10.0, "Ca9.0-s3p2",      "Ca9.0-s3p2d1",    "Ca9.0-s3p2d2",   ],
  ["Sc", "Sc_PBE19",  11.0, "Sc9.0-s2p2d1",    "Sc9.0-s3p2d1",    "Sc9.0-s3p2d2",   ],
  ["Ti", "Ti_PBE19",  12.0, "Ti7.0-s2p2d1",    "Ti7.0-s3p2d1",    "Ti7.0-s3p2d2f1", ],
  ["V",  "V_PBE19",   13.0, "V6.0-s2p2d1",     "V6.0-s3p2d1",     "V6.0-s3p2d2f1",  ],
  ["Cr", "Cr_PBE19",  14.0, "Cr6.0-s2p2d1",    "Cr6.0-s3p2d1",    "Cr6.0-s3p2d2f1", ],
  ["Mn", "Mn_PBE19",  15.0, "Mn6.0-s2p2d1",    "Mn6.0-s3p2d1",    "Mn6.0-s3p2d2f1", ],
  ["Fe", "Fe_PBE19H", 16.0, "Fe5.5H-s2p2d1",   "Fe5.5H-s3p2d1",   "Fe5.5H-s3p2d2f1",],
  ["Fe", "Fe_PBE19S", 14.0, "Fe6.0S-s2p2d1",   "Fe6.0S-s3p2d1",   "Fe6.0S-s3p2d2f1",],
  ["Co", "Co_PBE19H", 17.0, "Co6.0H-s2p2d1",   "Co6.0H-s3p2d1",   "Co6.0H-s3p2d2f1",],
  ["Co", "Co_PBE19S", 15.0, "Co6.0S-s2p2d1",   "Co6.0S-s3p2d1",   "Co6.0S-s3p2d2f1",],
  ["Ni", "Ni_PBE19H", 18.0, "Ni6.0H-s2p2d1",   "Ni6.0H-s3p2d1",   "Ni6.0H-s3p2d2f1",],
  ["Ni", "Ni_PBE19S", 16.0, "Ni6.0S-s2p2d1",   "Ni6.0S-s3p2d1",   "Ni6.0S-s3p2d2f1",],
  ["Cu", "Cu_PBE19H", 19.0, "Cu6.0H-s2p2d1",   "Cu6.0H-s3p2d1",   "Cu6.0H-s3p2d2f1",],
  ["Cu", "Cu_PBE19S", 11.0, "Cu6.0S-s2p1d1",   "Cu6.0S-s3p2d1",   "Cu6.0S-s3p2d2f1",],
  ["Zn", "Zn_PBE19H", 20.0, "Zn6.0H-s2p2d1",   "Zn6.0H-s3p2d1",   "Zn6.0H-s3p2d2f1",],
  ["Zn", "Zn_PBE19S", 12.0, "Zn6.0S-s2p1d1",   "Zn6.0S-s3p2d1",   "Zn6.0S-s3p2d2f1",],
  ["Ga", "Ga_PBE19",  13.0, "Ga7.0-s2p2d1",    "Ga7.0-s3p2d2",    "Ga7.0-s3p2d2f1", ],
  ["Ge", "Ge_PBE19",   4.0, "Ge7.0-s2p1d1",    "Ge7.0-s3p2d2",    "Ge7.0-s3p2d2f1", ],
  ["As", "As_PBE19",  15.0, "As7.0-s3p2d1",    "As7.0-s3p2d2",    "As7.0-s3p2d2f1", ],
  ["Se", "Se_PBE19",   6.0, "Se7.0-s3p2d1",    "Se7.0-s3p2d2",    "Se7.0-s3p2d2f1", ],
  ["Br", "Br_PBE19",   7.0, "Br7.0-s3p2d1",    "Br7.0-s3p2d2",    "Br7.0-s3p2d2f1", ],
  ["Kr", "Kr_PBE19",   8.0, "Kr10.0-s2p2d1",   "Kr10.0-s3p2d2",   "Kr10.0-s3p2d2f1",],
  ["Rb", "Rb_PBE19",   9.0, "Rb11.0-s2p2d1",   "Rb11.0-s3p2d2",   "Rb11.0-s3p2d2f1",],
  ["Sr", "Sr_PBE19",  10.0, "Sr10.0-s2p2d1",   "Sr10.0-s3p2d2",   "Sr10.0-s3p3d2f1",],
  ["Y",  "Y_PBE19",   11.0, "Y10.0-s3p2d1",    "Y10.0-s3p2d2",    "Y10.0-s3p3d2f1", ],
  ["Zr", "Zr_PBE19",  12.0, "Zr7.0-s3p2d1",    "Zr7.0-s3p2d2",    "Zr7.0-s3p2d2f1", ],
  ["Nb", "Nb_PBE19",  13.0, "Nb7.0-s3p2d1",    "Nb7.0-s3p2d2",    "Nb7.0-s3p2d2f1", ],
  ["Mo", "Mo_PBE19",  14.0, "Mo7.0-s3p2d1",    "Mo7.0-s3p2d2",    "Mo7.0-s3p2d2f1", ],
  ["Tc", "Tc_PBE19",  15.0, "Tc7.0-s3p2d1",    "Tc7.0-s3p2d2",    "Tc7.0-s3p2d2f1", ],
  ["Ru", "Ru_PBE19",  14.0, "Ru7.0-s3p2d1",    "Ru7.0-s3p2d2",    "Ru7.0-s3p2d2f1", ],
  ["Rh", "Rh_PBE19",  15.0, "Rh7.0-s3p2d1",    "Rh7.0-s3p2d2",    "Rh7.0-s3p2d2f1", ],
  ["Pd", "Pd_PBE19",  16.0, "Pd7.0-s3p2d1",    "Pd7.0-s3p2d2",    "Pd7.0-s3p2d2f1", ],
  ["Ag", "Ag_PBE19",  17.0, "Ag7.0-s3p2d1",    "Ag7.0-s3p2d2",    "Ag7.0-s3p2d2f1", ],
  ["Cd", "Cd_PBE19",  12.0, "Cd7.0-s3p2d1",    "Cd7.0-s3p2d2",    "Cd7.0-s3p2d2f1", ],
  ["In", "In_PBE19",  13.0, "In7.0-s3p2d1",    "In7.0-s3p2d2",    "In7.0-s3p2d2f1", ],
  ["Sn", "Sn_PBE19",  14.0, "Sn7.0-s3p2d1",    "Sn7.0-s3p2d2",    "Sn7.0-s3p2d2f1", ],
  ["Sb", "Sb_PBE19",  15.0, "Sb7.0-s3p2d1",    "Sb7.0-s3p2d2",    "Sb7.0-s3p2d2f1", ],
  ["Te", "Te_PBE19",  16.0, "Te7.0-s3p2d2",    "Te7.0-s3p2d2f1",  "Te7.0-s3p3d2f1", ],
  ["I",  "I_PBE19",    7.0, "I7.0-s3p2d2",     "I7.0-s3p2d2f1",   "I7.0-s3p3d2f1",  ],
  ["Xe", "Xe_PBE19",   8.0, "Xe11.0-s3p2d1",   "Xe11.0-s3p2d2",   "Xe11.0-s3p2d2f1",],
  ["Cs", "Cs_PBE19",   9.0, "Cs12.0-s3p2d1",   "Cs12.0-s3p2d2",   "Cs12.0-s3p2d2f1",],
  ["Ba", "Ba_PBE19",  10.0, "Ba10.0-s3p2d1",   "Ba10.0-s3p2d2",   "Ba10.0-s3p2d2f1",],
  ["La", "La_PBE19",  11.0, "La8.0-s3p2d1f1",  "La8.0-s3p2d2f1",  "La8.0-s3p3d2f1", ],
  ["Ce", "Ce_PBE19",  12.0, "Ce8.0-s3p2d1f1",  "Ce8.0-s3p2d2f1",  "Ce8.0-s3p3d2f1", ],
  ["Pr", "Pr_PBE19",  13.0, "Pr8.0-s3p2d1f1",  "Pr8.0-s3p2d2f1",  "Pr8.0-s3p3d2f1", ],
  ["Nd", "Nd_PBE19",  14.0, "Nd8.0-s3p2d1f1",  "Nd8.0-s3p2d2f1",  "Nd8.0-s3p3d2f1", ],
  ["Pm", "Pm_PBE19",  15.0, "Pm8.0-s3p2d1f1",  "Pm8.0-s3p2d2f1",  "Pm8.0-s3p3d2f1", ],
  ["Sm", "Sm_PBE19",  16.0, "Sm8.0-s3p2d1f1",  "Sm8.0-s3p2d2f1",  "Sm8.0-s3p3d2f1", ],
  ["Dy", "Dy_PBE19",  20.0, "Dy8.0-s3p2d1f1",  "Dy8.0-s3p2d2f1",  "Dy8.0-s3p3d2f1", ],
  ["Ho", "Ho_PBE19",  21.0, "Ho8.0-s3p2d1f1",  "Ho8.0-s3p2d2f1",  "Ho8.0-s3p3d2f1", ],
  ["Lu", "Lu_PBE19",  11.0, "Lu8.0-s3p2d2",    "Lu8.0-s3p2d2f1",  "Lu8.0-s3p3d2f1", ],
  ["Hf", "Hf_PBE19",  12.0, "Hf9.0-s3p2d2",    "Hf9.0-s3p2d2f1",  "Hf9.0-s3p3d2f1", ],
  ["Ta", "Ta_PBE19",  13.0, "Ta7.0-s3p2d2",    "Ta7.0-s3p2d2f1",  "Ta7.0-s3p3d2f1", ],
  ["W",  "W_PBE19",   12.0, "W7.0-s3p2d2",     "W7.0-s3p2d2f1",   "W7.0-s3p3d2f1",  ],
  ["Re", "Re_PBE19",  15.0, "Re7.0-s3p2d2",    "Re7.0-s3p2d2f1",  "Re7.0-s3p3d2f1", ],
  ["Os", "Os_PBE19",  14.0, "Os7.0-s3p2d2",    "Os7.0-s3p2d2f1",  "Os7.0-s3p3d2f1", ],
  ["Ir", "Ir_PBE19",  15.0, "Ir7.0-s3p2d2",    "Ir7.0-s3p2d2f1",  "Ir7.0-s3p3d2f1", ],
  ["Pt", "Pt_PBE19",  16.0, "Pt7.0-s3p2d2",    "Pt7.0-s3p2d2f1",  "Pt7.0-s3p3d2f1", ],
  ["Au", "Au_PBE19",  17.0, "Au7.0-s3p2d2",    "Au7.0-s3p2d2f1",  "Au7.0-s3p3d2f1", ],
  ["Hg", "Hg_PBE19",  18.0, "Hg8.0-s3p2d2",    "Hg8.0-s3p2d2f1",  "Hg8.0-s3p3d2f1", ],
  ["Tl", "Tl_PBE19",  19.0, "Tl8.0-s3p2d2",    "Tl8.0-s3p2d2f1",  "Tl8.0-s3p3d2f1", ],
  ["Pb", "Pb_PBE19",  14.0, "Pb8.0-s3p2d2",    "Pb8.0-s3p2d2f1",  "Pb8.0-s3p3d2f1", ],
  ["Bi", "Bi_PBE19",  15.0, "Bi8.0-s3p2d2",    "Bi8.0-s3p2d2f1",  "Bi8.0-s3p3d2f111",],
]
