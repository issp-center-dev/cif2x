.. _sec-cif2x-appendix:

================================================================
拡張ガイド
================================================================

mode を追加する
----------------------------------------------------------------

Quantum ESPRESSO の計算モードへの対応を追加するには、 ``src/cif2x/qe/calc_mode.py`` の ``create_modeproc()`` 関数に mode と変換クラスの対応付けを記述します。

.. code-block:: python

  def create_modeproc(mode, qe):
      if mode in ["scf", "nscf"]:
          modeproc = QEmode_pw(qe)
      else:
          modeproc = QEmode_generic(qe)
      return modeproc
	

modeごとの変換機能は ``QEmode_base`` の派生クラスとしてまとめられています。
このクラスは
``update_namelist()`` で namelist ブロックの更新と、
``update_cards()`` で cards ブロックのデータ生成を行います。
現在は pw.x の scf および nscf に対応する ``QEmode_pw`` クラスと、変換せずそのまま出力する ``QEmode_generic`` クラスが用意されています。

.. code-block:: python

  class QEmode_base:
      def __init__(self, qe):
      def update_namelist(self, content):
      def update_cards(self, content):


namelist については、空欄の値を結晶構造データ等から生成して代入するほか、格子パラメータなど Structure から決まる値や、他のパラメータとの整合性をとる必要のある値を強制的にセットする場合があります。処理内容はモードごとに個別に対応します。

cards ブロックについては、card の種類ごとに関数を用意し、card名と関数の対応付けを ``card_table`` 変数に列挙します。
基底クラスの ``update_cards()`` では、card名から対応する関数を取得して実行し、card の情報を更新します。もちろん、全く独自に ``update_cards()`` 関数を作成することもできます。

.. code-block:: python

    self.card_table = {
        'CELL_PARAMETERS': generate_cell_parameters,
        'ATOMIC_SPECIES': generate_atomic_species,
        'ATOMIC_POSITIONS': generate_atomic_positions,
        'K_POINTS': generate_k_points,
    }

cardごとの関数は ``scr/cif2x/qe/cards.py`` にまとめられており、関数名は ``generate_{card名}`` としています。この関数は card ブロックのパラメータを引数に取り、card名、option、dataフィールドからなる辞書データを返します。
