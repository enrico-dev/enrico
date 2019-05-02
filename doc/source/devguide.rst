Developer's Guide
=================

Contributing
------------

Any change to ENRICO must be made through a pull request. This applies to all
changes to documentation, code, binary files, etc. Even long term committers and
project owners must use pull requests.

No pull request may be merged without being independently reviewed.

Style Guidelines
----------------

.. important:: To ensure consistent styling with little effort, this project
    uses `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_. The
    repository contains a ``.clang-format`` file that can be used to
    automatically apply the style rules that are described below. The easiest
    way to use clang-format is through a plugin/extension for your editor/IDE
    that automatically runs clang-format using the ``.clang-format`` file
    whenever a file is saved.

Follow the `C++ Core Guidelines
<http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`_ except when they
conflict with another guideline listed here.

Conform to the C++14 standard.

Indentation
~~~~~~~~~~~

Use two spaces per indentation level.

Source Files
~~~~~~~~~~~~

Use a ``.cpp`` suffix for code files and ``.h`` for header files.

Header files should always use include guards with the following style (See
`SF.8 <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#sf8-use-include-guards-for-all-h-files>`_):

.. code-block:: C++

    #ifndef ENRICO_MODULE_NAME_H
    #define ENRICO_MODULE_NAME_H

    namespace enrico {
    ...
    content
    ...
    }

    #endif // ENRICO_MODULE_NAME_H

Avoid hidden dependencies by always including a related header file first,
followed by local includes, other library includes, and then C/C++ library
includes. For example:

.. code-block:: C++

   // foo.cpp
   #include "foo.h"

   #include "bar.h"
   #include "baz.h"

   #include <gsl/gsl>

   #include <cstddef>
   #include <iostream>
   #include <vector>

Naming
~~~~~~

Struct and class names should be CamelCase, e.g. ``HexLattice``.

Functions (including member functions) should be lower-case with underscores,
e.g. ``get_indices``.

Local variables, global variables, and struct/class member variables should be
lower-case with underscores, except for physics symbols that are written
differently by convention (e.g., ``E`` for energy). Data members of classes
additionally have trailing underscores (e.g., ``a_class_member_``).

All classes and non-member functions should be declared within the ``enrico``
namespace.

Accessors and mutators (get and set functions) may be named like
variables. These often correspond to actual member variables, but this is not
required. For example, ``int count()`` and ``void set_count(int count)``.

Variables declared constexpr or const that have static storage duration (exist
for the duration of the program) should be upper-case with underscores,
e.g., ``SQRT_PI``.

Use C++-style declarator layout (see `NL.18
<http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#nl18-use-c-style-declarator-layout>`_):
pointer and reference operators in declarations should be placed adjacent to the
base type rather than the variable name. Avoid declaring multiple names in a
single declaration to avoid confusion:

.. code-block:: C++

   T* p; // good
   T& p; // good
   T *p; // bad
   T* p, q; // misleading

Source Annotation
~~~~~~~~~~~~~~~~~

Classes, structs, and functions are to be annotated for the `Doxygen
<http://www.stack.nl/~dimitri/doxygen/>`_ documentation generation tool. Use the
``\`` form of Doxygen commands, e.g., ``\brief`` instead of ``@brief``.

Doxygen
-------

The Doxygen documentation (including a full description of the class hierarchies) can be found
`here <doxygen/html/index.html>`_.
