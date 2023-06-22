#ifndef _TEST_H_
#define _TEST_H_

#include "ParticlePositions.h"

namespace unfolding {
	class test {
	public:
		ParticlePositions _partpos;
		test(ParticlePositions PartPos): _partpos(PartPos){};
		void printnumber() {
			cout << _partpos.n_part << endl;
		}
	};
}


#endif
