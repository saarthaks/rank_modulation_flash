# Rank Modulation Codes for NAND Flash Memory

In this project, I've explored the potential for a relatively new type of ECC for storing data in NAND Flash memory. Several papers since about 2008 have been developing these codes, which use permutations over $[n]$ to encode information rather than the $q$-ary codes we've seen in class. From a coding perspective, the main idea is to mitigate challenges in quantizing the underlying analog state of a Flash cell and instead using that analog information for more precise comparative orderings over the charge-state of a Flash cell. This work explores some theoretical guarantees that the current literature has proven on the existence and bounds of their rate and distance. It then returns to the motivating domain of error-correcting Flash memory and explores some Monte-Carlo simulations of these error-correcting codes on a simplified model of Flash cell retention loss.

# References
[1] Barg, A., and Mazumdar, A. "Codes in permutations and error correction for rank mod-ulation." _IEEE Inter. Symp. on Info. Theory_ (2010).

[2] Cai, Yu, et al. “Data retention in MLC NAND Flash memory: Characterization, opti-mization, and recovery.” _IEEE Inter. Symp. on High Perf. Computer Architecture_ (2015).

[3] Cai, Yu, et al. “Error analysis and retention-aware error management for NAND Flash memory.” _Intel Technology Journal_ (2013).

[4] Cai, Yu, et al. "Threshold voltage distribution in MLC NAND flash memory: Characterization, analysis, and modeling." _Design, Automation, and Test in Europe Conference Exhibition (DATE)_ 2013.

[5] Chadwick, H., and Kurz, L. "Rank permutation group codes based on Kendall’s correla-tion statistic." _IEEE Trans. on Info. Theory_ 15.2 (1969): 306-315.

[6] Farnoud, F., Skachek, V., and Milenkovic, O. "Error-correction in flash memories via codes in the Ulam metric." _IEEE Trans. on Info. Theory_ 59.5 (2013): 3003-3020.

[7] Gabrys, R. et. al. "Codes correcting erasures and deletions for rank modulation." _IEEE Trans. on Info. Theory_ 62.1 (2015): 136-150.

[8] Jiang, A., Schwartz, M., and Bruck, J. "Error-correcting codes for rank modulation." _IEEE Inter. Symp. on Info. Theory_ (2008).

[9] Luo, Yixin, et al. "Improving 3D NAND flash memory lifetime by tolerating early retention loss and process variation." _Proceedings of the ACM on Measurement and Analysis of Computing Systems_ (2018).

[10] Mazumdar, A., Barg, A., and Zemor, G. "Constructions of rank modulation codes." _IEEE Trans. on Info. Theory_ 59.2 (2012): 1018-1029.

[11] Yehezkeally, Y., and Schwartz, M. "Snake-in-the-box codes for rank modulation." _IEEE Trans. on Info. Theory_ 58.8 (2012): 5471-5483.
