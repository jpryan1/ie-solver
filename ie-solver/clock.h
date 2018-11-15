#ifndef _CLOCK_H_
#define _CLOCK_H_

namespace ie_solver{

struct Clock{
	double time;
	double elapsed_;

	void tic(){
		time = omp_get_wtime();
	}
	void toc(){
		elapsed_ += (omp_get_wtime()-time);
	}
	void elapsed(const char* s){
		printf("%s: %f seconds elapsed\n", s, elapsed_);
	}
	void toc(const char* s){
		printf("%s: %f seconds\n", s, (omp_get_wtime() - time));
	}
};

} // namespace

#endif