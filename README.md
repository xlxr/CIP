[![INFORMS Journal on Computing](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# [On exact approaches for interdicting connected dominating sets](https://doi.org/10.1287/ijoc.2025.1541)

This archive accompanies the paper *On exact approaches for interdicting connected dominating sets* by Haonan Zhong, Foad Mahdavi Pajouh, and Oleg A. Prokopyev in the INFORMS Journal on Computing. It contains the C++ implementations of the two main algorithms and the 135 weighted graph instances. This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).
## Cite

Please cite both the paper and this software archive using their respective DOIs.

Paper — <https://doi.org/10.1287/ijoc.2025.1541>

Software and data — <https://doi.org/10.1287/ijoc.2025.1541.cd>

BibTeX for this snapshot of the repository follows.

```bibtex
@misc{ZhongEtAl2026CIP,
  author    = {Zhong, Haonan and Mahdavi Pajouh, Foad and Prokopyev, Oleg A.},
  publisher = {INFORMS Journal on Computing},
  title     = {{On exact approaches for interdicting connected dominating sets}},
  year      = {2026},
  doi       = {10.1287/ijoc.2025.1541.cd},
  url       = {https://github.com/INFORMSJoC/2025.1541},
  note      = {Available for download at https://github.com/INFORMSJoC/2025.1541}
}
```

## Description

The connected dominating set interdiction problem asks for a minimum-cost set of edges to block so that every connected dominating set in the remaining graph has vertex weight at least the threshold `eta`.

| Source file | Algorithm | Paper | Output |
| --- | --- | --- | --- |
| `src/IP.cpp` | Lazy branch-and-cut (LBC) | Section 3.2 | `IP.csv` |
| `src/CBB.cpp` | Combinatorial branch-and-bound (CBB) | Section 4 | `CBB.csv` |

## Building

The experiments used 64-bit Windows, an Intel Core i9-14900HX, 16 GB of RAM, one core, and a 3,600-second limit per instance. CBB uses the C++ standard library. LBC also requires the Gurobi C++ library and a Gurobi license. The paper used Gurobi Optimizer 10.0.3.

From the repository root in an **x64 Native Tools Command Prompt for Visual Studio**, build CBB with

```bat
cl /nologo /O2 /EHsc /std:c++17 src\CBB.cpp /FeCBB.exe
```

With Gurobi 10.0.3 installed and `GUROBI_HOME` pointing to its `win64` directory, build LBC in the same prompt. This example uses the release C++ library with the dynamic Visual C++ runtime. [Gurobi 10 Windows installation guide](https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/10.0/quickstart_windows.pdf)

```bat
cl /nologo /O2 /EHsc /std:c++17 /MD /I"%GUROBI_HOME%\include" src\IP.cpp /FeIP.exe /link /LIBPATH:"%GUROBI_HOME%\lib" gurobi_c++md2017.lib gurobi100.lib
```

## Replicating

Each program reads `filename.txt` in its working directory. The first line is the number of instances. Each remaining line names one graph file. The following Windows commands run one instance from Table 3 in a separate directory.

```bat
mkdir run-ieee30
copy data\real\IEEE-30-0.333333.graph run-ieee30\
cd run-ieee30
>filename.txt echo 1
>>filename.txt echo IEEE-30-0.333333.graph
..\CBB.exe
..\IP.exe
```

For this instance, Table 3 reports `eta = 55`, an optimal interdiction cost of **2** for both algorithms, and a minimum-cut cost of **8**. In `CBB.csv`, the corresponding columns are `Best Sol`, `Best Bound`, and `mincut`, with values 2, 2, and 8.

| Results in the paper | Directory containing the graphs and `filename.txt` | Instances |
| --- | --- | ---: |
| Table 1, random graphs with 30, 50, or 70 vertices | `data/random_small/` | 45 |
| Table 2, random graphs with 100, 200, or 300 vertices | `data/random_large/` | 45 |
| Table 3, real-life networks | `data/real/` | 45 |

For a complete table, open a prompt at the repository root and switch to the corresponding directory. For example, run Table 1 with

```bat
cd data\random_small
..\..\CBB.exe
..\..\IP.exe
```

Each executable writes a fresh CSV in its working directory when it starts. Save that CSV before repeating the same run. When running LBC on Linux, use LF line endings in `filename.txt` so its entries match the graph filenames. Graph formats and sources are documented in the [data README](data/README.MD).

## Development and support

The authors maintain the [CIP development repository](https://github.com/xlxr/CIP). Contact information is in [AUTHORS](AUTHORS).
