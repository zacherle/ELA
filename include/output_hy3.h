#ifndef OUTPUT_HY3_H
#define OUTPUT_HY3_H

void output_hy3(struct hy3_file &hy3, TParams &param, const double hypo[4], 
		const CArrivalFile &arrs, std::vector<TRecordHyp> &hyp, Gather &gather,
		const double* co, int info);

#endif // OUTPUT_HY3_H
