.. SWGOMP documentation master file, created by
   sphinx-quickstart on Wed Apr 17 13:59:04 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##########
SWGOMP手册
##########

什么是SWGOMP？
--------------

SWGOMP是基于swgcc的GOMP开发的OpenMP Offload功能，包括编译器插件和运行时库，可以支持Fortran、C、C++语言。

接下来应该？
-------------

如果需要快速使用SWGOMP，请参考 :ref:`quickstart`。

了解SWGOMP的工作原理，请参考 :ref:`how_worked`。

使用SWGOMP遇到了问题，可参见 :ref:`faq`。

声明
----

许可证
~~~~~~

SWGOMP的插件源码使用 `GPLv3 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_ 授权许可，如果基于插件部分的代码进行改动后进行发布，则应遵守GPL授权协议（主要是继续开源与专利相关事宜）。

SWGOMP的运行时库源码使用 `LGPLv3 <https://www.gnu.org/licenses/lgpl-3.0.en.html>`_ 授权，链接SWGOMP运行时库并进行发布，则不必继续开源。

**注意** ，两种协议都声明了::

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

或者说，本项目 **不提供任何代码质量保证** 。

使用SWGOMP完成的科研论文，建议引用该文献： `Lin, H., Yan, L., Chang, Q. et al. O2ath: an OpenMP offloading toolkit for the sunway heterogeneous manycore platform. CCF Trans. HPC (2024). <https://doi.org/10.1007/s42514-024-00191-1 >`_ 。

开发者
~~~~~~

段晓辉（山东大学，sunrise.duan@sdu.edu.cn）

####
目录
####

.. toctree::
   :maxdepth: 2

   quickstart
   how_worked
   faq
