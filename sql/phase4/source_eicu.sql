-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
BEGIN READ ONLY;
\copy (SELECT cellattributepath,celllabel,cellattribute,count(*) AS records FROM public.nursecare WHERE concat_ws(' ',cellattributepath,celllabel,cellattribute) ~* 'arteri|radial|femoral|brachial|axillar|catheter|line' GROUP BY 1,2,3) TO '@CUFF_WORKSPACE@/phase4/site_nursecare_dictionary_eicu.csv' CSV HEADER
\copy (SELECT cplgroup,cplitemvalue,count(*) AS records FROM public.careplangeneral WHERE concat_ws(' ',cplgroup,cplitemvalue) ~* 'arteri|radial|femoral|brachial|axillar|catheter|line' GROUP BY 1,2) TO '@CUFF_WORKSPACE@/phase4/site_careplangeneral_dictionary_eicu.csv' CSV HEADER
ROLLBACK;
