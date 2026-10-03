function [w,wp,wpp,wppp] = eval_aifantis_reference_1d(ref,r)
%EVAL_AIFANTIS_REFERENCE_1D Evaluate the 1-D reference and radial derivatives.
%
% Values outside [RI,RE] are clamped to the reference interval.  This is
% convenient later when straight-sided tetrahedral boundary facets slightly
% approximate the circular cylinder.

r = min(max(r,ref.r(1)),ref.r(end));
w    = ref.wInterp(r);
wp   = ref.wpInterp(r);
wpp  = ref.wppInterp(r);
wppp = ref.wpppInterp(r);
end
