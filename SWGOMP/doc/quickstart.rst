.. _quickstart:

快速开始
########

**注意** ，本项目目前仅兼容静态链接，使用 ``-cgsp 64 -b`` 提交的应用。

本节简要介绍如何开始使用SWGOMP，但由于神威平台的复杂性，一些细节在本节中被隐藏，了解这些细节显然有助于你在使用过程中定位错误，如果需要了解请移步 :ref:`how_worked`。也有一些常见的问题可以在 :ref:`faq` 中找到答案。

**建议通读本章后再上手**

获取源码
--------

源码位于： `<https://gitee.com/swmore/swgomp>`_

完成SWGOMP的编译
----------------

获取源码后直接在源码顶层目录执行：

.. code-block:: bash

  make

编译完成将生成运行时库 ``lib/libswgomp.a`` 和插件 ``plugin/swgomp.so`` 。

**注意**，SWGOMP将会使用宿主机的 ``g++`` 编译插件，并使用 ``sw9gcc`` 编译运行时库，请确认你的环境中这两个命令是存在的。

在构建工具中替换你的编译器
--------------------------

由于种种限制，我最后只能采用编译器包装脚本的方式提供带有OpenMP的编译一条龙服务，这些脚本位于 ``script`` 目录，其中：

===========    ===================
脚本           用途
===========    ===================
xcc.py         单节点C
xcxx.py        单节点C++
xfort.py       单节点Fortran
mpixcc.py      MPI C
mpixcxx.py     MPI C++
mpixfort.py    MPI Fortran
===========    ===================

你应当在构建工具中替换对应的神威平台编译器，从而让你的OpenMP代码可以被从核并行执行。
你可以使用绝对路径进行替换，也可以在

.. code-block:: bash

  source <你的SWGOMP>/script/setenv
  
之后直接使用脚本的名称。

另外，SWGOMP与SWACC一个比较大的区别是，SWGOMP使用了GCC的前端，所以不需要 ``.rmod`` 文件，对于Fortran库可以直接使用已存在的 ``.mod`` 文件而不需要重新使用SWGOMP编译库（需要在从核上调用该库除外）。故一开始你不需要重新编译你的Fortran库。


在你动手之前
------------

与普通的对称多核OpenMP不同，OpenMP Offload需要使用 ``target`` 语句将代码放到计算设备端，再使用 ``parallel`` 语句在计算设备端并行执行。

由于OpenMP Offload的设计中过多地考虑了GPU， ``target`` 提供了对应的 ``teams`` 语句，一个 ``team`` 类似于GPU上的一个 ``block``。 不同 ``team`` 内的线程老死相往来。

SWGOMP确实支持 ``teams distribute`` 那一套，但是感觉 **确实没有什么用** 。

本文中的例子将以C和Fortran为主，因为C++的话可以使用 `SWUC <https://arxiv.org/abs/2208.00607>`_ 套件，那样更加优雅一些。而且SWGOMP对虚函数没有特别好的支持

修改你的代码
------------

第一个并行循环
~~~~~~~~~~~~~~

当你并不知道SWGOMP如何工作时，你只要知道一个简单的从核并行循环如下所示：

Fortran：

.. code-block:: Fortran

  subroutine vfma(a, b, c, n)
    implicit none
    real(kind=8), intent(inout) :: a(n)
    real(kind=8), intent(in)    :: b(n), c(n)
    integer(kind=4), intent(in) :: n
    integer(kind=4)             :: i

    !$omp target parallel do
    do i = 1, n
      a(i) = a(i) + b(i)*c(i)
    end do
    !$omp end target parallel do
  end subroutine vfma

C：

.. code-block:: C

  void vfma(double *a, double *b, double *c, int n) {
    #pragma omp target parallel for
    for (int i = 0; i < n; i++) {
      a[i] += b[i]*c[i];
    }
  }

其实OpenMP Offload与OpenMP最大的区别就是，Offload多了一层 ``target`` 。

Fortran也有不需要循环的数组整体操作语法，可以使用OpenMP的 ``workshare`` 语句进行并行，例如：

.. code-block:: Fortran

  subroutine vfma(a, b, c, n)
    implicit none
    real(kind=8), intent(inout) :: a(n)
    real(kind=8), intent(in)    :: b(n), c(n)
    integer(kind=4), intent(in) :: n
    integer(kind=4)             :: i

    !$omp target
    !$omp parallel workshare
    a(:) = a(:) + b(:)*c(:)
    !$omp end parallel workshare
    !$omp end target
  end subroutine vfma

但是gfortran前端会给使用OpenMP语句带来一些限制，详见 :ref:`ffe_bug` 一节。

线程私有变量
~~~~~~~~~~~~

有一些循环使用了外部变量，这些变量并不真的会产生依赖，例如：

.. code-block:: C

  void vfma(double *a, double *b, double *c, int n) {
    double prod; //说的就是你
    #pragma omp target parallel for
    for (int i = 0; i < n; i++) {
      prod = b[i]*c[i];
      a[i] += prod;
    }
  }

但是语义上似乎每个线程都在向 ``prod`` 写入，从而产生不可预料的结果。
这里的 ``prod`` 就是一个需要改为线程似有的变量（或许也不需要，我有种印象是变量默认 ``private`` ，但是为了更完善的声明，不然你也处理一下？）

对于C语言，进行变量私有化的方式有两种，第一种是利用C语言特性，将变量移到内部的作用域：

.. code-block:: C

  void vfma(double *a, double *b, double *c, int n) {
    #pragma omp target parallel for
    for (int i = 0; i < n; i++) {
      double prod = b[i]*c[i]; //这样prod会开在线程的栈上
      a[i] += prod;
    }
  }

另一种是利用OpenMP的 ``private`` 子句：

.. code-block:: C

  void vfma(double *a, double *b, double *c, int n) {
    double prod;
    #pragma omp target parallel for private(prod) //声明prod为线程似有
    for (int i = 0; i < n; i++) {
      prod = b[i]*c[i];
      a[i] += prod;
    }
  }

Fortran则只能使用 ``private`` 子句，因为Fortran没有嵌套的作用域：

.. code-block:: Fortran

  subroutine vfma(a, b, c, n)
    implicit none
    real(kind=8), intent(inout) :: a(n)
    real(kind=8), intent(in)    :: b(n), c(n)
    integer(kind=4), intent(in) :: n
    integer(kind=4)             :: i
    real(kind=8)                :: prod !Fortran的定义只能在函数开头，很讨厌
  
    !$omp target parallel do private(prod) //使用private语句的例子
    do i = 1, n
      prod = b(i)*c(i)
      a(i) = a(i) + prod
    end do
    !$omp end target parallel do
  end subroutine vfma

注意，这里将变量声明为 ``private`` 的子句其实不会将变量的初始值拷入私有变量中，如果需要变量的初始值，可以使用 ``firstprivate`` 子句。

特别需要注意的是，应当私有化的变量漏掉的时候如果运气好并不会产生太坏的影响（例如编译器帮你加了private），但是 **数组一般不会被自动私有化** ，这是需要重点检查的情况。
一般来说，我通过检查Fortran的变量定义中，数组有没有与循环下标对应的维度来确定是否需要进行私有化。

输出的变量和需要归约的变量
~~~~~~~~~~~~~~~~~~~~~~~~~~

**注意：目前没有研究GOMP归约数组的机制，也没有对应地实现，但我明确知道归约数组和归约变量在GOMP中是两套逻辑**

归约变量可以使用 ``parallel`` 子句的 ``reduction`` 子句，team之间不能进行归约。同时，也应注意到，归约的结果是写到team的首从核的，所以还需要进行 ``map`` 操作。例如：

.. code-block:: Fortran

  !$omp target map(tofrom: sum) !启动时将主核的sum拷贝到副本，结束时将副本拷回主核的sum
  !$omp parallel do reduction(+:sum) !归约到0号从核的sum副本
  do ie = 1, n
    sum = sum + a(i)
  end do
  !$omp end parallel do
  !$omp end target

**注意** ，由于浮点加法不满足 **结合律** ，所以使用 ``reduction`` 会引入误差。

**注意** ，其他的并行区输出变量的情况也要用 ``tofrom`` 。

处理跨文件调用
~~~~~~~~~~~~~~

SWGOMP具有自动为被Offload的代码段生成从核版本的功能，但是仅限于编译当前源文件时可见的函数，如果Offload调用了其他源文件的函数，则SWGOMP无法进行自动化的处理。（呼叫从核LTO功能）

这种情况需要手动为被调用的函数添加从核化标记，具体语法为：

C语言：

.. code-block:: c
  
  int f(int a){
    return a;
  }
  //第一种方式，为单个函数声明从核化
  #pragma omp declare target(f)
  //第二种方式，为一段代码中的函数都声明从核化
  #pragma omp declare target
  int g(int a){
    return a;
  }
  int h(int a){
    return a;
  }
  #pragma omp end declare target

Fortran语言：

.. code-block:: Fortran

  subroutine f(a)
    implicit none
    !$omp declare target !大约可以加在这里，应该在定义变量的代码那里都可以
    real(kind=8), intent(inout) :: a

    return a
  end subroutine

.. _ffe_bug:

Fortran前端的一些限制
~~~~~~~~~~~~~~~~~~~~~

注意，gfortran-7.1.0的前端有一些bug，故你不能写：

.. code-block:: Fortran

  !$omp target parallel
  !$omp do
  !!你的循环在这里
  !$omp end do
  !$omp end target parallel

会导致报 ``Internal Compiler Error`` 错误， 幸运的是，你仍然可以写：

.. code-block:: Fortran

  !$omp target
  !$omp parallel do
  !!你的数组操作在这里
  !$omp end parallel do
  !$omp end target

同样， ``!$omp target parallel workshare`` 也是不被允许的，但是你只要将 ``target`` 单独拆出来就能正常使用了。
这件事情非常蠢，但是无法在编译器插件的层面解决。

.. _move_your_stack:

迁移你的栈
----------

默认情况下，私有变量位于LDM中，当你需要大量的私有变量时，LDM很有可能不够用。
此时，常见的报错为从核报错 ``Accessed Wrong Address`` ，并且使用 ``addr2line`` 定位PC会发现报错的行在函数的定义附近。
必要时可以为 ``target`` 语句添加 ``device(1)`` 子句，则可以将从核栈迁移到一个默认每从核8MiB的从核共享栈：

.. code-block:: Fortran

  subroutine vfma(a, b, c, n)
    implicit none
    real(kind=8), intent(inout) :: a(n)
    real(kind=8), intent(in)    :: b(n), c(n)
    integer(kind=4), intent(in) :: n
    integer(kind=4)             :: i
    real(kind=8)                :: prod(131072) !Fortran的定义只能在函数开头，很讨厌
  
    !$omp target parallel do private(prod) device(1) //不加device的话LDM会不够用哦
    do i = 1, n
      prod(1) = b(i)*c(i)
      a(i) = a(i) + prod
    end do
    !$omp end target parallel do
  end subroutine vfma

如果你的程序需要大量内存，无法为从核提供每从核8MiB的共享栈，或者你的并行区栈非常深，8MiB仍不能满足从核需求，可以通过设置 ``OMP_STACKSIZE`` 环境变量为从核提供你希望的贡献栈大小。例如：

.. code-block:: bash

  OMP_STACKSIZE=8M bsub -b -cgsp 64 ...
  OMP_STACKSIZE=512K bsub -b -cgsp 64 ...

此处的 ``M`` 和 ``K`` 都是2基的字节数，或者说， ``1M=1024*1024`` ， ``1K=1024`` 。

