# lammps文档

不用windows了，直接用linux版本

要得到自组装的图片即可

## 零 软件使用

### 1.如何打开示例

使用lammps gui和lammps

lammps gui打开方式

```
cd ~/chem/LAMMPS_GUI
./lammps-gui &
```

lammps gui打开例程

Flie-open打开in.xxxx文件

点击lammps-gui左下角run的按钮即可运行并查看结果

![](/home/rqxx/.config/marktext/images/2026-02-03-01-39-25-image.png)

运行结果有各类数据的随时间变化曲线

在右上角可以更改数据类型

Temp温度

<u>E_pair对势能：**系统中所有非键结原子对之间的相互作用势能在，自组装过程中，疏水尾部会聚集以减少与溶剂的接触，这通常会导致 `E_pair`显著下降并趋于稳定。如果这个值没有明显降低，可能意味着自组装没有发生或模拟时间不够。</u>**

E_mod分子势能

TolEng总能量

Press压力

### 2.可视化

使用Ovito可视化软件

```
cd /home/rqxx/chem/ovito-basic-3.14.1-x86_64/bin
./ovito 
```

打开可视化软件后

file->Load file打开run后生成的dump.xxxx文件即可查看到分子轨迹

### 3.输入文件如何编写

```
# 2d micelle simulation

dimension    2

neighbor    0.3 bin
neigh_modify    delay 5

atom_style    bond

# Soft potential push-off

read_data    data.micelle
special_bonds    fene

pair_style    soft 1.12246
pair_coeff    * * 0.0 1.12246

bond_style     harmonic
bond_coeff    1 50.0 0.75

velocity    all create 0.45 2349852

variable    prefactor equal ramp(1.0,20.0)

fix        1 all nve
fix        2 all temp/rescale 100 0.45 0.45 0.02 1.0
fix        3 all adapt 1 pair soft a * * v_prefactor
fix        4 all enforce2d

thermo        50
run        1000

unfix        3

# Main run

pair_style    lj/cut 2.5

# solvent/head - full-size and long-range

pair_coeff    1 1 1.0 1.0 2.5
pair_coeff    2 2 1.0 1.0 2.5
pair_coeff    1 2 1.0 1.0 2.5

# tail/tail - size-averaged and long-range

pair_coeff    3 3 1.0 0.75 2.5
pair_coeff    4 4 1.0 0.50 2.5
pair_coeff    3 4 1.0 0.67 2.5

# solvent/tail - full-size and repulsive

pair_coeff    1 3 1.0 1.0 1.12246
pair_coeff    1 4 1.0 1.0 1.12246

# head/tail - size-averaged and repulsive

pair_coeff    2 3 1.0 0.88 1.12246
pair_coeff    2 4 1.0 0.75 1.12246

thermo        50

#dump        1 all atom 2000 dump.micelle

#dump        2 all image 2000 image.*.jpg type type zoom 1.6
#dump_modify    2 pad 5 adiam 1 0.5 adiam 2 1.5 adiam 3 1.0 adiam 4 0.75

#dump        3 all movie 2000 movie.mpg type type zoom 1.6
#dump_modify    3 pad 5 adiam 1 0.5 adiam 2 1.5 adiam 3 1.0 adiam 4 0.75

reset_timestep    0
run        1000
```

示例代码如上

参数说明：



输入文件需要定义分子之间的健信息

### 4.data文件

data文件一般用于定义原子的位置

构建data文件的方式为

采用目前已有的可视化软件（例如material studio)先建模，然后将建好的模型，采用Lammps的软件模块msi2lmp转化为lammps格式的data文件

具体流程如下：

MS模型建好后，选择Moduals->Forcite->Calculation，选择cvff⼒场

计算结果输出选择.car格式，⽣成.car和.mdf⽂件

*data文件的建模是需要把分子都放入一个盒子然后进行的*

把.car .mdf  文件放到，mis2lmp下新建的model文件夹内

msi2lmp [⽂件名前缀] -class I -frc cvff -i

*（ -i 是为了忽略缺失化学键等信息的报错直接生成.data文件，注意这里我们只需要吧原子坐标的信息弄下来，所以可以忽略这种报错）*

生成.data文件（生成文件中的一些错误可以手动修改）

这样生成的data文件会有健信息的缺失



根据官方文档，如果为不过于复杂体系的data文件，可以自己用xxxx方法进行建立



## 一 方案选择

1.粗粒化模拟

2.全原子模拟
