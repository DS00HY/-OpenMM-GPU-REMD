# REMD

本程序是基于OpenMM实现的跨节点GPU并行REMD框架。

## 目录

- [使用说明](#使用说明)
- [待开发问题](#待开发问题)


### 使用说明

#### 1.安装
在使用之前确保当前环境已经包含下面工具
- MPICH: 版本高于3.3
- SWIG
- Anaconda，并且安装好Python的Numpy、MPI4PY、OpenMM包
- GCC
- Make
- CUDA

#### 2.C++部分编译

- 在源码文件输入命令“make build”，便开始执行编译过程。
- 如果没有提示错误，目录下生成mpiplus_wrap.cxx, mpiplus.py,_mpiplus.so三个文件则说明编译成功，
  c++的扩展代码可以当做Python包正常使用了。
- 如果想删除编译的代码，可以输入“make clean”。
- 在Makefile的第22行，可以设置是否开启debug模式，输出C++的信息，如果不想开启debug模式，将22行的${DEBUG}注释掉即可。
- 如果担心C++代码的编译情况，可以写一个简单的测试脚本来测试C++部分是否可以正常运行，
  可以以附录中test.py为例进行测试。将test.py代码运行，如果结果为没有运行错误，
  并且查看C++和Python代码输出的paramlist和id相互对应，即可说明C++代码编译成功，并且功能正常。
 （其中test.py其26行中getExchange_func()函数是mpiplus模块提供给Python的调用接口，comm为MPI.COMM, i表示当前的迭代次数，state表示当前进程的需要传递给主进程的状态数组，paramlist为进程的参数列表，ex_kind和md_kind为后续开发选项，目前值只支持(1,1)一种。） 

#### 3.REMD使用

- 以4个副本的testREMD.py的脚本跨节点运行的资源分配命令如下：
   python setconfigfile.py --np 4 python testREMD.py -debug
   其中-np 指定运行的副本数目，testREMD.py为需要分配节点资源运行的REMD程序，-debug命令为是否输出相关的debug信息。该命令完成后，会生成configfile 和hostfile两个文件来保存进程的运行信息。
- 将代码运行到指定节点和GPU上，在脚本中继续输入以下命令即可：
     mpiexec.hydra -f hostfile -configfile configfile
     实际运行的脚本可以参考run.sh/ rerun.sh的写法
- 结果输出：以样例testREMD.py的测试结果的输出为例:
   PDB	out/out0.pdb, out/out1.pdb.... 
   模拟信息	out/md0.log, out/md1.log....
   Debug信息	debug.log
   每次交换对象	ex_pair.out
   温度副本参数	replica.out
   断点文件	Checkpoint0.44 (副本0，第44次迭代的断点)
   其中值得说明的是，replica.out的最后一行为下一次迭代的次序，还有每个副本应该使用的温度列表的id，因此该文件也是断点重启的重要依据。
#### 4.断点重启
如果需要断点重启，在原来程序的基础上，将REMD.run()的restart参数设置为True即可开始继续运行。需要注意的是其模拟结果的输出文件最好进行编号，以副本0为例，进行第一次重启的pbd可以命名为out0-1.pdb，第二次为out0-2.pdb。如果采用这种命名方式，可以使用提供的mergepdb.py进行pdb的整合。将所有的out0.pdb, out0-1.pdb, ... , out0-n.pdb以及保存重启信息的restart.log文件与mergepdb程序放到一起，使用mergerpdb中定义的函数get_pdb(pdb_numb, rank, pdb_iteration, n_iteration)，即可将多次重启的pdb文件进行整合，该函数的参数pdb_numb表示需要整合的pdb数目，rank表示当前副本号，pdb_iteration为运行代码中的pdb存储间隔，n_iteration表示运行的总次数。
实际运行的程序可以参考restartREMD.py, 提交脚本为rerun.sh的写法。


### 待开发问题

- 让用户指定温度连续还是轨迹连续
- 结果pdb的压缩存储
- 与Gromacs匹配

