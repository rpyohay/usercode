% gg_to_gg.mp -- do not edit, generated automatically by my_gluon.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 40*2.84526pt#, 25*2.84526pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=2.84526pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1, __i2);
vright (__o1, __o2);
vlabel ("$g$", __i1);
vlabel ("$g$", __i2);
vlabel ("$g$", __o1);
vlabel ("$g$", __o2);
vconnect ("curly", __i1, __v1, __i2);
vconnect ("curly", __o1, __v2, __o2);
vconnect ("curly,label=$g$,label.dist=4mm", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
