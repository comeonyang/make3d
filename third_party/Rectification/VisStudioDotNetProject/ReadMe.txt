========================================================================
    CONSOLE APPLICATION : VisStudioDotNetProject Project Overview
========================================================================

See USAGE in the directory above for details on the code.

It's the most commonly asked question so I'll state for the record that you should take care when linking to your version of lapack.
It's necessary to make sure the runtime libraries are the same for both the clapack.lib & your code.
Depends on how you've compiled clapack - but mine is on Multi-thredaded Debug DLL so you may need to change it.

Both debug & release targets work so take your pick.
