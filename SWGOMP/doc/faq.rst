.. _faq:

常见问题
########

使用
====

我有一个非常非常长的 ``!$omp`` 制导语句，如何对其进行换行？
-----------------------------------------------------------

对于 free form 的 Fortran 代码，这里可以采用的换行方式主要有两种，第一种是使用 ``&`` 续行，续行后继续 ``!$omp`` ， 例如：

.. code-block:: Fortran

  !$omp parallel private(a, b, c, d) &
  !$omp private(e, f, g, h)
  a = omp_get_thread_num()
  !$omp end parallel

第二种是在开启了 ``-cpp`` 选项的情况下，使用 ``\`` 进行续行。

.. code-block:: Fortran

  !$omp parallel private(a, b, c, d \
                         e, f, g, h)
  a = omp_get_thread_num()
  !$omp end parallel

如何为主从核采用不同的编译选项？
--------------------------------

``swgomp-driver.py`` 中提供了 ``-Whost`` 和 ``-Wslave`` 选项，可以单独为主核或从核提供传给SWGCC的选项，例如：

.. code-block:: bash

  sw9gcc -Whost,-fno-inline -Wslave,-funroll-loops -c test.c

这样可以使主核编译时加入 ``-fno-inline`` ，在从核编译时加入 ``-funroll-loops`` 。

报错
====

Accessed Wrong Address/Unaligned Exception
------------------------------------------

请使用 ``addr2line`` 定位，如果在函数定义的行，一种常见的情况是LDM栈溢出，可以参考 :ref:`move_your_stack` 。

Internal Compiler Error
-----------------------

很抱歉会引发这种问题，如果遇到这种问题，可以通过以下步骤确认问题是源自SWGOMP还是GOMP还是原生编译器：

1. 使用原生编译器替代编译包装脚本来编译你的程序，如果可以通过，那么说明问题来自GOMP或者SWGOMP。
2. 使用原生编译器 + ``-fopenmp`` 替代编译包装脚本来编译你的程序，如果可以通过，说明问题来自SWGOMP，不要犹豫提交一个issue。否则，略微拆分你的OpenMP源语试一下，可能是 ``gfortran`` 的GOMP支持部分有问题，详见 :ref:`ffe_bug` 。

性能
====

部分循环在从核上非常非常慢
--------------------------

经我的调研，该问题主要来自于LDCache颠簸。一般来说，这些循环可能涉及到了非常大量的数组，如果恰巧有 ``>4`` 个数组映射到了同一个Cache块，那么LDCache的命中率是神仙难救的。

据此，我提出两种解决方案：

第一种针对 ``>4`` ，也就是采用 ``omnicopy`` 取一部分数据，降低使用LDCache的数组个数。

第二种针对同一个Cache块，可以在分配数组地址时使它们的地址在 ``mod(32768)`` 的意义下均匀分布。这要求DIY一个内存分配器。我针对GRIST写过一个，但目前没有严苛测试，有需要者可以与我联系。

兼容性
======

我是否还能使用athread库？
-------------------------

你只可以使用普通的 ``athread_spawn`` 和 ``athread_join`` 。并且仅在 ``cgsp 64`` 的情况下使用。

我是否还能使用LDM和DMA？
------------------------

并行区默认的栈会放在LDM。

提供了 ``omnicopy`` 函数以支持DMA。C语言中该函数接口与 ``memcpy`` 相同，但是会根据源地址和目标地址自动切换为DMA方式取数据。

Fortran语言的 ``omnicopy`` 函数由 ``omnicopy.mod`` 提供，支持整型和浮点型数组。接口为：

.. code-block:: Fortran

  call omnicopy(target, src) !拷贝与target大小相当的数据
  call omnicopy(target, src, count) !count个数据，这里不需要乘数据类型大小了
  call omnicopy_c(target, src, size) !拷贝大小为size字节

我可以在动态链接或者akernel程序中使用SWGOMP吗？
-----------------------------------------------

目前来看不行，因为SWGOMP使用了wrap对athread库进行了包装，这些在动态链接中无法使用，如果有迫切需求请提交issue。

过去基于xfort开发的程序能使用SWGOMP吗？
---------------------------------------

xfort太丑了，并且考虑到其并没有公开的release，所以在设计中没有考虑xfort的兼容性。如果希望将xfort的程序迁移到SWGOMP，可以参考的改动如下：

.. code-block:: Fortran

  !$omp target
  call get_coreid(tid)
  do i = tid, n, 64
    !...
  end do
  !$omp end target

迁移到SWGOMP只需为其再添加 ``parallel`` 语句，然后将 ``get_coreid`` 替换为标准的 ``omp_get_thread_num`` （其中需要引入 ``omp_lib`` 模块）。

.. code-block:: Fortran

  use omp_lib !引用omp_lib获得omp_get_thread_num函数

  !$omp target
  !$omp parallel !添加parallel语句
  tid = omp_get_thread_num()
  do i = tid, n, 64
    !...
  end do
  !$omp end parallel
  !$omp end target

