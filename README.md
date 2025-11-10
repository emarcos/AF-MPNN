# AFM-MPNN and AF3-MPNN: AlphaFold predictions guided by designed MSAs

Accurate prediction of synthetic protein–protein interactions is critical for designing novel therapeutics, diagnostics, and imaging reagents. Deep learning has revolutionized protein structure prediction, yet modeling heteromeric and loop-mediated interactions—as in antibodies, nanobodies, and monobodies—remains challenging due to limited evolutionary information at interfaces and the conformational flexibility of loops. Here, we show that AlphaFold-Multimer (AFM) and AlphaFold3 (AF3) can be guided using paired multiple sequence alignments (MSAs) generated from ProteinMPNN interface design. AFM and AF3 with these designed MSAs (AFM-MPNN and AF3-MPNN) focus sampling around the input model, overcoming default limitations. Benchmarking on datasets of nanobody– and monobody–receptor complexes demonstrates that designed MSAs encode interface contacts and guide predictions toward accurate binding orientations and CDR loop conformations. Improvements are specific to binding complexes, with lower accuracy and confidence for non-binding decoys, highlighting their potential for discriminating binders from non-binders. AF3-MPNN also enables structural refinement of partially correct models and improves decoy ranking in challenging, low-confidence de novo predictions. Its computational efficiency and direct compatibility with MSA-driven deep learning tools make it broadly applicable for both protein design and structure prediction workflows, opening avenues for the de novo design of complex and irregular binding interfaces previously inaccessible.

This computational approach, its rationale and benchmarks across multiple datasets of nanobodies and monobodies are presented in the article: "Designed MSAs guide AlphaFold for improved modeling of loop-mediated protein-protein interactions" by Lourdes Carcelén, Alexandre Casadesús and Enrique Marcos.

Dependencies:
- ColabFold
- ProteinMPNN
- PyRosetta
- TM-align
