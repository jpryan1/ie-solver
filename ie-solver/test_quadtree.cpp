#include <stdlib.h>
#include <stdio.h>
#include "quadtree.h"

int main(int argc, char** argv){

	printf("Hello world\n");
	QuadTree q;

	std::vector<float> arr(50);
	for(int i=0; i<50; i++){
		arr[i] = (i+1)/51.0;
	}
	q.initialize_tree(arr);
	q.print();
	return 0;
}