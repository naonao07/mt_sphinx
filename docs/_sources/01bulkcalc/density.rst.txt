密度の算出方法
=====================
液体の密度はNPT計算を行うことで算出出来ます。  


まずは水の密度を求めてみます。

水分子の用意
--------------------


はじめに、計算する水分子のパラメータを設定します。

今回はSPC/E力場パラメータを使って密度を計算してみましょう。

水の力場とそれを使ったlammpsスクリプトがLAMMPSのマニュアルに記載があるので、この内容を参考にして書いていきます。

`SPC water model (lammps ドキュメント) <https://docs.lammps.org/Howto_spc.html>`_


::

    units real
    atom_style full
    region box block -5 5 -5 5 -5 5
    create_box 2 box  bond/types 1 angle/types 1 &
                    extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2

    mass 1 15.9994
    mass 2 1.008

    pair_style lj/cut/coul/cut 10.0
    pair_coeff 1 1 0.1553 3.166
    pair_coeff 1 2 0.0    1.0
    pair_coeff 2 2 0.0    1.0

    bond_style zero
    bond_coeff 1 1.0

    angle_style zero
    angle_coeff 1 109.47


1行目はスクリプト内で使う単位系を設定しています。

今回使っているrealという単位系は、例えば質量の単位が :math:`g/mol` なので、7行目のmassコマンドで設定している数値は :math:`15.9996 g/mol` を意味しています。 
`units (lammps ドキュメント) <https://docs.lammps.org/units.html>`_



次に2行目ですが、これは原子一個一個にどんな値を設定できるか(atom_style)を決めています。このスクリプトでのatom_styleはfullなので、次のような数値を設定することが可能です。設定しなくても構いません。

    * atom-ID ：          各原子の名前
    * type    ：          識別ID
    * position ：         座標 x,y,z
    * velosity ：         速度 vx,vy,vz
    * force    ：         働く力 fx,fy,fz
    * image-frags：       可視化のための数値データ
    * mask       ：       所属するグループ名
    * bond        ：      2原子間の力場パラメータデータ
    * angle        ：     3原子間の力場パラメータデータ
    * dihedral     ：     4原子間の力場パラメータデータ
    * improper     ：     以上以外の構造における力場パラメータデータ
    * charge      ：      各原子の電荷

`atom_style (lammps ドキュメント) <https://docs.lammps.org/atom_style.html>`_

このページでは特別な理由がない限り*unit real*を用いることにします。



次に3行目、この行ではregionコマンドを使って計算するセルボックスの大きさを設定しています。この設定ではちょうど以下の画像で描いた領域が計算セルボックスになります。

.. image:: fig01.png

コードの意味は、左から 

region 領域名 領域の型 型に合わせた設定値

となります。

`region (lammps ドキュメント) <https://docs.lammps.org/region.html>`_ を見ると





次の4行目は、上で用意した計算セルにどれくらいの原子種(atom_type)、結合種(bond_type)、・・・を入れることにするかを設定しています。