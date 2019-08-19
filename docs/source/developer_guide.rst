Developer Guide
===============

``stk`` is a fairly simple library. The key design principle is
modularity and atomicity. The basic dependency structure is that
any module can depend on things defined in the :mod:`~.stk.utilities`
module but all other top-level modules, such as :mod:`.calculators`
and :mod:`.molecular` should be completely independent of one another.

Within these top level modules, sub-modules should still be defined in
as dependency free way as possible, though exceptions can be made in
rare and extreme circumstances where it makes sense. For example
:mod:`.molecules` imports a couple of tools from
:mod:`.functional_groups`


:mod:`~.stk.utilities` defines general purpose tools, :mod:`.molecular`
defines tools dealing with molecules and :mod:`.calculators` defines
calculator objects. Each calculator type should be defined in its
own sub-module with additional sub-folders to make everything as tidy
and logically cohesive as possible.
