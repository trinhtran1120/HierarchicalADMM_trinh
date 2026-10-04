---
title: Notes for EPCN ProjectHybrid ADMM

---

# Notes for EPCN Project
## Hybrid ADMM
### 1. Previous interpretation
- Previously, we considered hybrid ADMM as a hierarchical combination of two ADMM variants: one algorithm is applied at the global level, while another is applied within local modules
![image](https://hackmd.io/_uploads/rkLruOltMx.png)
***Figure 1. Illustration of 2 types of hybrid ADMM algorithms***. 

- Each level executes its own ADMM procedure, and we mainly intervene in how these levels interact.

### 2. Dr. Binh's idea
- A more flexible definition is to construct the hybrid algorithm by intervening the iterations and communication queries of fADMM and hADMM.

- Instead of assigning one ADMM variant to an entire level of the hierarchy, we can selectively use fADMM or hADMM updates for different parts of the tree. 
For example, module 2 can use hADMM while module 3 uses fADMM. This provides freedom for restructuring the computational tree and, in particular, for making tree more balanced.

Consider the following tree

![hybrid_ADMM-Page-1.drawio](https://hackmd.io/_uploads/H1QDYuxFMx.png)
***Figure 2. New representation of hybrid ADMM***.


- Instead of grouping Modules 2 and 3 into a single computational module like Section 1 (Figure 1),we can separate them into distinct algorithmic modules, as in Figure 2.

- Note that if Module 2 requires information from Module 3, Module 3 must still be queried before the corresponding update of Module 2 can be completed.

![hybrid_ADMM-Page-2.drawio](https://hackmd.io/_uploads/B1r2K_xFfx.png)
***Figure 3. Intuition behind the proposed hybrid ADMM***.

<mark>**The key difference from the previous approach is that hybridization now occurs within the iteration/query structure, rather than only between nested inner and outer ADMM algorithms.**</mark>