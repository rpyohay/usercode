% gg_to_tt_to_WbWb_to_ffbenub_to_ffbenuba.mp -- do not edit, generated automatically by my_gluon.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 80*2.84526pt#, 50*2.84526pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=2.84526pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1, __i2);
vright (__o1, __o2, __o3, __o4, __o5, __o6, __o7);
vlabel ("$g$", __i1);
vlabel ("$g$", __i2);
vlabel ("$\bar {f}^{'}$", __o1);
vlabel ("$f$", __o2);
vlabel ("$b$", __o3);
vlabel ("$\bar {b}$", __o4);
vlabel ("$e^{+}$", __o5);
vlabel ("$\gamma $", __o6);
vlabel ("$\nu _{e}$", __o7);
vconnect ("curly", __i1, __v1, __i2);
vconnect ("fermion", __o3, __v3);
vconnect ("fermion,label=$\bar {t}$,label.side=left", __v3, __v4);
vconnect ("fermion,label=$t$", __v4, __v5);
vconnect ("fermion", __v5, __o4);
vconnect ("curly,label=$g$,label.dist=4mm", __v1, __v4);
vconnect ("photon,label=$W^{+}$,label.side=left", __v5, __v7);
vconnect ("photon,label=$W^{-}$", __v3, __v2);
vconnect ("fermion", __o1, __v2, __o2);
vconnect ("fermion", __o5, __v6);
vconnect ("fermion", __v6, __v7);
vconnect ("photon,tension=0", __v6, __o6);
vconnect ("fermion", __v7, __o7);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
