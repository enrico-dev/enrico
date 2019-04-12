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
followed by C/C++ library includes, other library includes, and then local
includes. For example:

.. code-block:: C++

   // foo.cpp
   #include "foo.h"

   #include <cstddef>
   #include <iostream>
   #include <vector>

   #include "hdf5.h"
   #include "pugixml.hpp"

   #include "error.h"
   #include "random_lcg.h"

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
pointer and reference operators in declarations should be placed adject to the
base type rather than the variable name. Avoid declaring multiple names in a
single declaration to avoid confusion:

.. code-block:: C++

   T* p; // good
   T& p; // good
   T *p; // bad
   T* p, q; // misleading

Curly braces
~~~~~~~~~~~~

For a class declaration, the opening brace should be on the same line that
lists the name of the class.

.. code-block:: C++

    class Matrix {
      ...
    };

For a function definition, the opening and closing braces should each be on
their own lines.  This helps distinguish function code from the argument list.
If the entire function fits on one or two lines, then the braces can be on the
same line. e.g.:

.. code-block:: C++

    return_type function(type1 arg1, type2 arg2)
    {
      content();
    }

    return_type
    function_with_many_args(type1 arg1, type2 arg2, type3 arg3,
                            type4 arg4)
    {
      content();
    }

    int return_one() {return 1;}

    int return_one()
    {return 1;}

For a conditional, the opening brace should be on the same line as the end of
the conditional statement. If there is a following ``else if`` or ``else``
statement, the closing brace should be on the same line as that following
statement. Otherwise, the closing brace should be on its own line. A one-line
conditional can have the closing brace on the same line or it can omit the
braces entirely e.g.:

.. code-block:: C++

    if (condition) {
      content();
    }

    if (condition1) {
      content();
    } else if (condition 2) {
      more_content();
    } else {
      further_content();
    }

    if (condition) {content()};

    if (condition) content();

For loops similarly have an opening brace on the same line as the statement and
a closing brace on its own line. One-line loops may have the closing brace on
the same line or omit the braces entirely.

.. code-block:: C++

    for (int i = 0; i < 5; i++) {
      content();
    }

    for (int i = 0; i < 5; i++) {content();}

    for (int i = 0; i < 5; i++) content();

Documentation
~~~~~~~~~~~~~

Classes, structs, and functions are to be annotated for the `Doxygen
<http://www.stack.nl/~dimitri/doxygen/>`_ documentation generation tool. Use the
``\`` form of Doxygen commands, e.g., ``\brief`` instead of ``@brief``.
