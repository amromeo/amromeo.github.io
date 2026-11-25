# Task G decision

Accepted 25 November 2026. Method A passed all 16 prespecified conditions on the 7,488 independent Stage 3 confirmation datasets. It had lower all-dataset loss than the checked blend at each N and met every margin for low-pathology loss and paired error. The fixed-runtime audit reproduced all 749 fits. [Full results](results/REPORT.md)

On descriptive RIbench at ≤20% pathology, the blend-first order had slightly lower loss at N = 1,000, 5,000 and 50,000. The paired 95% intervals included zero at all three sizes. The package proposal keeps the prespecified order at every N. Run safeguarded anchoredRI first. If it declines or withholds, try the Stage 3 post-fit blend (t = 0.25 on the Stage 2 peak-choice fit) and release its limits only if the Stage 1 check passes. [Protocol](PROTOCOL.md)
