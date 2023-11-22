# Personalized Pricing with Group Fairness Constraint

## Introduction
This repository explores a personalized pricing strategy under fairness constraints, using total variation distance (TVD) and earth mover's distance (EMD) as fairness measures. We introduce a linear programming approach in discrete price settings and provide theoretical bounds for the revenue gap between discrete and continuous price settings. Our model demonstrates enhanced fairness in terms of customer surplus, social welfare, and Gini index while maintaining competitive pricing policies.

## Repository Structure
- `code/`: Contains all the R scripts implementing analyzing the proposed and baseline models.
- `result/`: Includes results such as optimal objective, solution, and evaluation metrics for linear demand models.
- `results/`: Houses results for other demand models including optimal objectives, solutions, and evaluation metrics.

## Files Overview
- `GroupFair_EM_FlowLP_VD_*`: Scripts related to models with fairness constraints measured by EMD.
- `GroupFair_TV_FlowLP_VD_*`: Scripts associated with models constrained by TVD.
- `GroupFair_constr_meanp_FlowLP_VD_*`: Scripts for models with fairness and mean price constraints.
- `GroupFair_unconstr_FlowLP_VD_*`: Unconstrained models for comparison.

## Usage
Each script is named according to the model and the constraints it implements. Run the scripts in the `code/` directory to reproduce the results presented in the paper.

## Contributing
Contributions are welcome! If you have any suggestions or improvements, please feel free to fork the repository and submit a pull request.

## License
Please refer to the LICENSE file for details on the use and distribution of this repository's contents.

## Citation
If you use the code or data from this repository in your research, please cite our paper as follows:
```bibtex
@inproceedings{10.1145/3593013.3594097,
author = {Chen, Xin and Xu, Zexing and Zhao, Zishuo and Zhou, Yuan},
title = {Personalized Pricing with Group Fairness Constraint},
year = {2023},
isbn = {9798400701924},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3593013.3594097},
doi = {10.1145/3593013.3594097},
booktitle = {Proceedings of the 2023 ACM Conference on Fairness, Accountability, and Transparency},
pages = {1520–1530},
numpages = {11},
keywords = {personalized pricing, social welfare, group fairness, statistical parity},
location = {Chicago, IL, USA},
series = {FAccT '23}
}
```
## Contact
For any questions or further discussion regarding the research, please open an issue on this repository.

