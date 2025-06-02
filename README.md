# Fill‑in‑the‑Middle Code‑Completion Benchmark

**Author:** Daniel Dajbov

**Project:** Personal C++ Huffman encoder project, evaluated with *StarCoder2‑3B*

**Goal:** Build a small, *explainable* test bed (≈40 prompts) that measures how well modern code‑completion models can fill in the gap a programmer is about to write.

---

## 1  Project Overview

This repo walks through the **whole evaluation loop**:

1. **Dataset creation** – slice a real‑world C++ file (`src/huffman_code.cpp`) into 40 *prefix / middle / suffix* triples that mimic the user’s cursor position.
2. **Automatic scoring** – run any Hugging Face model with a single function call and compute two complementary metrics:
   * **TSED** – *Tree‑Sitter Edit Distance* (syntax‑aware similarity, 0 – 1).
   * **chrF** – character n‑gram F‑score (surface similarity, 0 – 100).
3. **Manual sanity check** – inspect each prediction/ground‑truth pair and decide which metric agrees better with human judgement.

Everything is bundled into **`evaluation_pipeline.ipynb`** – open the notebook, hit **Run‑All**, and you’ll finish with plots, CSV‑ready scores, and ground‑truth vs. model snippets saved in `predictions_ground_truth/`.

---

## 2  Repository Layout

```text
Fill-in-the-middle-problem-main/
├── src/
│   └── huffman_code.cpp            # C++ file that seeds the dataset
├── grammar.json                    # C++ Tree‑Sitter grammar (pruned)
├── evaluation_pipeline.ipynb       # dataset ▸ inference ▸ metrics ▸ plots
├── predictions_ground_truth/
│   ├── ground_truth_sample_XX.txt  # canonical middle snippets
│   └── predicted_code_sample_XX.txt# model completions
└── README.md                     
```

## 3  Running the Pipeline

1. Launch `jupyter lab`, open `evaluation_pipeline.ipynb`.
2. Change **`API_URL`** to test another model (any FIM‑capable causal LM).
3. **Run all cells**.


## 4  Metric Rationale

| Metric   | What it checks                                                                                                                          | Why it matters                                                                        |
| -------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| **TSED** | Normalised tree‑edit distance between ASTs (uses C++ Tree‑Sitter). Rewards syntactic correctness & penalises structural hallucinations. | Code that *parses/compiles* is usually the first requirement for a useful completion. |
| **chrF** | Character‑level F‑score (n‑gram overlap). Robust to tokenisation differences, very fast.                                                | Complements TSED by catching *semantic* tokens even when structure matches.           |

Manual inspection showed TSED correlates **better** with perceived correctness (≈ 0.72 Spearman vs. chrF’s ≈ 0.55).



## 5  Manual Validation Workflow

A lightweight rubric was used:

1. **Perfect (✓)** – model snippet identical in behaviour.
2. **Acceptable (\~)** – compiles & passes eye‑test but uses different naming / minor style.
3. **Wrong (✗)** – syntax error *or* wrong semantics.


## 6  License
Code released under the **MIT License**.
