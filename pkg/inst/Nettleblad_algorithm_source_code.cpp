// cnF2freq, (c) Carl Nettelblad, Department of Information Technology, Uppsala University 2008
// Public release 0.4
//
// carl.nettelblad@it.uu.se
//
// This code is allowed to be freely used for any commercial or research purpose. If the code is
// used or integrated into another project largely unchanged, attribution to the original author
// is appreciated. No warranties are given.


#define NDEBUG
// These defines fixed an error in one particular site installation of the Portland compiler.
#define _STLP_EXPOSE_GLOBALS_IMPLEMENTATION 1
#define _REENTRANT 1
#define _SECURE_SCL 0
// For MSCVC
// Note: MSVC OpenMP support is insufficient in current release. Disable OpenMP for compilation
// in MSVC.
//
// Recent releases of g++ on MacOS X and Linux, as well as the Intel C++ compiler on
// Linux (x86 and x64) and Windows have been tested. The Portland compiler collection works only
// with some optimization settings, some race conditions in STL, although some
// workarounds are used.
#define _CRT_SECURE_NO_WARNINGS

#include <vector>

#include <stdio.h>


// _MSC_VER is here to be interpreted as any compiler providing TR1 C++ headers
#ifdef _MSC_VER
#include <array>
#else
// Boost also provides an array implementation, that is largely compatible
#include <boost/array.hpp>
#endif

#include <stdlib.h>
#include <set>
#include <algorithm>
#include <math.h>
#include <float.h> // use these libraries


using namespace std; // use functions that are part of the standard library
#ifdef _MSC_VER
using namespace tr1;
#else
using namespace boost;
#endif

#ifndef _MSC_VER
#define _isnan isnan
#define _finite finite
#endif


const int sexmarkerval = 66;// make sure this doesn't occur in the data
// F2 with haplotyping
/*const int NUMGEN = 3;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS] = {0, 0, 1, 1, 0, 1};
const int NUMTYPES = 1 << TYPEBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = NUMGEN;
const unsigned int HALFNUMPATHS = 1 << (TYPEBITS / 2);
const unsigned int NUMPATHS = NUMTYPES << 1;
const unsigned int NUMSHIFTGEN = NUMGEN - 1;
const unsigned int HALFNUMSHIFTS = 1 << ((1 << (NUMSHIFTGEN - 1)) - 1);
const unsigned int NUMSHIFTS = 1 << ((1 << NUMSHIFTGEN) - 1); 
const bool HAPLOTYPING = true;*/


// F2 with no haplotyping
const int NUMGEN = 2;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS] = {0, 1};
const int NUMTYPES = 1 << TYPEBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = 1;
const unsigned int HALFNUMPATHS = 1;
const unsigned int NUMPATHS = 2;
const unsigned int NUMSHIFTGEN = 0;
const unsigned int HALFNUMSHIFTS = 1;
const unsigned int NUMSHIFTS = 1;
const bool HAPLOTYPING = false;



const int HALFNUMTYPES = 1 << (TYPEBITS / 2);

// Infer corrections for impossible genotypes according to the pedigree
// Also infer values for missing markers from existing information in offspring
// When enabled, results similar to QTL Express without data reduction
// ccoeff does not provide correction inference, so exact result reproduction
// is achieved when this flag is disabled.
const bool CORRECTIONINFERENCE = true;

bool early = false;

vector<double> markerposes;
vector<double> actrec[2];
vector<unsigned int> chromstarts;

vector<int> markertranslation;

double discstep = 0.1;
double baserec[2];
int sexc = 1;


// Is a specific marker value a admissible as a match to marker b
// The logic we use currently represents "0" as unknown, anything else as a known
// marker value.
// NOTE, VALUE CHANGED FOR ZEROS!
template<bool zeropropagate> bool markermiss(unsigned int& a, const unsigned int b)
{
	// This is the real logic; we do not bind anything at all when zeropropagate is true
	if (zeropropagate) return false;

	if (!a)
	{
		if (!zeropropagate) a = b;
		return false;
	}
	if (!b && a != sexmarkerval) return false;

	return a != b;
}


// Move a tree-based binary flag up a generation. The structure of bit flags might look like
// (1)(3)(3), where each group of (3) is in itself (1)(2), where (2) is of course (1)(1).
int upflagit(int flag, int parnum, int genwidth)
{
	if (flag < 0) return flag;
	flag >>= parnum * (genwidth - 1);
	flag &= ((1 << (genwidth - 1)) - 1);

	return flag;
}

struct individ;

int generation = 1;
int shiftflagmode;
int lockpos[NUMSHIFTS];
int quickmark[NUMSHIFTS];
int quickgen[NUMSHIFTS];

// The quick prefixes are caches that retain the last invocation.
// Only used fully when we stop at marker positions exactly, i.e. not a general grid search.
double quickfactor[NUMSHIFTS];
array<double, NUMTYPES> quickendfactor[NUMSHIFTS];
array<array<double, NUMTYPES>, NUMTYPES> quickendprobs[NUMSHIFTS];
array<double, NUMTYPES> quickmem[NUMSHIFTS];

// A hashed store of inheritance pathway branches that are known to be impossible.
// Since we can track the two branches that make up the state in the F_2 individual independently,
// this optimization can reduce part of the cost by sqrt(number of states).
typedef array<array<array<array<array<array<array<int, 4>, HALFNUMSHIFTS>, HALFNUMPATHS + 1>, HALFNUMTYPES>, 2>, 2>, 2> IAT;
IAT impossible;

// A memory structure storing haplo information for later update.
// By keeping essentially thread-independent copies, no critical sections have to
// be acquired during the updates.
array<array<float, 2>, 2000000> haplos;

// done, factors and cacheprobs all keep track of the same data
// done indicates that a specific index (in the binary tree of blocks of multi-step transitions) is done
// with a "generation id" that's semi-unique, meaning no active clearing of the data structure is performed
vector<int> done[NUMSHIFTS];
// factors contain the mantissas of the extended floating-point representation
vector<array<float, NUMTYPES> > factors[NUMSHIFTS];
// cacheprobs contain actual transitions from every possible state to every possible other state
vector<array<array<float, NUMTYPES>, NUMTYPES> > cacheprobs[NUMSHIFTS];
vector<individ*> reltree;

//#pragma omp threadprivate(realdone, realfactors, realcacheprobs)


#pragma omp threadprivate(generation, done, factors, cacheprobs, shiftflagmode, impossible, haplos, lockpos, quickmark, quickgen, quickmem, quickfactor, quickendfactor, quickendprobs, reltree)

// We put all thread-local structures in our own separate struct. This is because many compilers implementing OpenMP
// use a relatively expensive call for determining thread-local global data, which is not cached over function calls.
// The simple pointer-arithmetic for a lookup within a struct is cheaper.
struct threadblock
{
	int* const generation;
	int* const shiftflagmode;
	int* const quickmark;
	int* const quickgen;
	int* const lockpos;
	double* const quickfactor;
	array<double, NUMTYPES>* const quickmem;
	IAT* const impossible;
	array<array<float, 2>, 2000000>* const haplos;
	vector<int>* const done;
	vector<array<float, NUMTYPES> >* const factors;
	vector<array<array<float, NUMTYPES>, NUMTYPES> >* const cacheprobs;
	array<double, NUMTYPES>* const quickendfactor;
	array<array<double, NUMTYPES>, NUMTYPES>* const quickendprobs;

	threadblock() : generation(&::generation), shiftflagmode(&::shiftflagmode), impossible(&::impossible),
		done(::done), factors(::factors), cacheprobs(::cacheprobs), haplos(&::haplos),
		quickmark(::quickmark), quickgen(::quickgen), lockpos(::lockpos), quickmem(::quickmem),
		quickfactor(::quickfactor), quickendfactor(::quickendfactor), quickendprobs(::quickendprobs)
	{
	};
};


// Turners are used as mix-ins when probabilities are filtered at the "fixed" marker
// (the one we are actually probing).
// This functionality is used to invert the state in different manner, to test for the
// possibility that the complete haplotype assignment from an arbitrary point until the end
// has been mixed up. This can easily happen as the assignment of haplotype numbers is
// basically arbitrary, so the only thing keeping it in place is the linkage to the original
// defining locus.
class noneturner
{
public:
	void operator() (const double* probs) const
	{
	};
} none;

class aroundturner
{
	int turn;
	int flagmodeshift;

public:
	aroundturner(int turn) : turn(turn & 54), flagmodeshift(
		(turn >> TYPEBITS)
		| ((turn & 1) ? 2: 0)
		| ((turn & 8) ? 4 : 0))
	{
	}


	void operator() (double* probs) const
	{
		double probs2[NUMTYPES];
		for (int i = 0; i < NUMTYPES; i++)
		{
			probs2[i ^ turn] = probs[i];
		}

		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] = probs2[i];
		}

		// Not *tb out of laziness, the number of calls is limited
		shiftflagmode ^= flagmodeshift;
	};
};

// A struct containing some auxiliary arguments to the trackpossible family of functions.
// These parameters are not known at compile-time, hence not template arguments, but they do not
// change over the series of recursive calls.
const struct trackpossibleparams
{
	const float updateval;
	int* const gstr;

	trackpossibleparams() : updateval(0.0f), gstr(0)
	{
	}

	trackpossibleparams(float updateval, int* gstr) : updateval(updateval), gstr(gstr)
	{
	}
} tpdefault;


// A structure containing most of the information on an individual
struct individ
{
	// The individual #.
	int n;
	// Generation number. Convenient, while not strictly needed.
	int gen;
	// Parents.
	individ* pars[2];
	// Sex.
	bool sex;
	// Line or strain of origin, should only exist in founders.
	int strain;
	// Marker data as a list of pairs. No specific ordering assumed.
	vector<pair<int, int> > markerdata;
	// Temporary storage of all possible marker values, used in fixparents.
	vector<set<int> > markervals;
	// The haplotype weight, or skewness. Introducing an actual ordering of the value in markerdata.
	vector<float> haploweight;
	// The cost-benefit value of inverting the haplotype assignment from an arbitrary marker point on.
	vector<float> negshift;

	// Accumulators for haplotype skewness updates.
	vector<float> haplobase;
	vector<float> haplocount;

	individ()
	{
		pars[0] = 0;
		pars[1] = 0;
		strain = 0;
		sex = false;
	}


	// A wrapper for calling trackpossible. Additional logic introduced to do lookup in the "impossible" tables of branches
	// not allowed at the current point. This will generally mean most branches and reduces the number of states visited considerably
	// in typical cases.
	//
	// The double overload means that the class can be treated as a conventional function call, returning a double, when the pre-lookup
	// is not needed. In the "real" trackpossible method, one instance is called in that way, while the other one is pre-looked up.
	template<bool update, bool zeropropagate> struct recursetrackpossible
	{
		individ* mother;
		int upflagr;
		int upflag2r;
		int upshiftr;
		const trackpossibleparams& extparams;
		const int genwidth;
		const int markerval;
		const int marker;
		const threadblock& tb;
		int firstpar;
		int* impossibleref;
		int impossibleval;
		bool prelok;

		recursetrackpossible(individ* mother, const threadblock& tb, int markerval, int marker, int upflag, int upflag2, int upshift, int genwidth, int f2n, int firstpar, int numrealpar, const trackpossibleparams& extparams) :
		genwidth(genwidth), extparams(extparams), markerval(markerval), marker(marker), tb(tb), firstpar(firstpar), mother(mother)

		{
			upflagr = upflagit(upflag, firstpar, genwidth);
			upflag2r = upflagit(upflag2, firstpar, genwidth >> (NUMGEN - NUMFLAG2GEN));
			upshiftr = upflagit(upshift, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));

			prelok = true;
			if (!zeropropagate && genwidth == (1 << (NUMGEN - 1)))
			{
				impossibleref = &(*tb.impossible)[*(tb.shiftflagmode) & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3];
				impossibleval = (*tb.generation) * markerposes.size() + marker;

				if (*impossibleref == impossibleval)
				{
					prelok = false;
				}
			}
		}

		operator double()
		{
		  if (!prelok)
		  {
			  return 0;
		  }

			double baseval =
				mother->pars[firstpar]->trackpossible<update, zeropropagate>(tb, markerval, marker,
				upflagr,
				upflag2r,
				upshiftr, extparams, genwidth >> 1);

			if (!zeropropagate && !update && genwidth == (1 << (NUMGEN - 1)) && !baseval)
			{
				*impossibleref = impossibleval;
			}

			return baseval;
		}
	};

	// zeropropagate also implies that there is a gstr value, we propagate zeros to find any possible source strain
	// The main logic of tracking a specific inheritance pathway, computing haplotype weights, and overall feasibility.
	// update: Should haplotype weights be updated?
	// zeropropagate: Should zero marker values be kept, or changed to the actual values they are matched against.
	// threadblock: Reference to all thread-local data.
	// inmarkerval: The marker val we are matching against.
	// marker: The marker number.
	// flag: The actual genotype flag. Note that this is zero-extended to be one bit more than the actual enumeration of flags
	// (lowest bit always 0).
	// flag99: What's mostly called flag2. A specific shift state of marker numbers, or -1 to indicate that all shifts are allowed.
	// localshift: Based on shiftflagmode, the overall mapping *for the complete sequence* of strand numbers to actual parentage for all
	// individuals in the analysis.
	// extparams: external parameters.
	// genwidth: The width of the generation flags.
	template<bool update, bool zeropropagate> double trackpossible(const threadblock& tb, unsigned int inmarkerval, const unsigned int marker,
		const unsigned int flag, const int flag99, int localshift = 0, const trackpossibleparams& extparams = tpdefault,
		const int genwidth = 1 << (NUMGEN - 1)) /*const*/
	{
		if (this == NULL) return 1;

		int upflag2 = -1;
		const int upflag = flag >> 1;
		const int upshift = localshift >> 1;
		int f2s = 0;
		const int* themarker = &markerdata[marker].first;
		bool allthesame = themarker[0] == themarker[1];
		int f2end = 2;

		if (flag99 != -1 && genwidth >> (NUMGEN - NUMFLAG2GEN) > 0)
		{
			upflag2 = flag99 >> 1;
			f2s = flag99;
			f2end = flag99 + 1;			
		}

		int firstpar = flag & 1;
		double ok = 0;

#pragma ivdep
		// This ivdep is quite complicated, since we actually change markerval, but not in
		// coniditions where this alters the results of the other iteration.

		// flag2 determines which value in the tuple to check. flag and localshift determine the haplotype weight value
		// assigned to that test, and which parent to match against that value.
		for (int flag2 = f2s; flag2 < f2end && (HAPLOTYPING || !ok); flag2++)
		{

			unsigned int markerval = inmarkerval;
			double baseval;

			int f2n = (flag2 & 1);
			int realf2n = f2n;

			// If this marker value is not compatible, there is no point in trying.
			if (markermiss<zeropropagate>(markerval, themarker[f2n])) continue;

			// Normalize, in some sense.
			f2n ^= ((firstpar ^ localshift) & 1);

			if (!genwidth)
			{
				//printf("Hoj!\n");
			}

			if (zeropropagate || !genwidth)
			{
				baseval = 0.5;
			}
			else if (allthesame)
			{
				baseval = (f2n) ? 1.0 : 0.0;
			}
			else
			{
				if (HAPLOTYPING)
					baseval = fabs((f2n ? 1.0 : 0.0) - haploweight[marker]);
				else
				{
					// No haplotype weights, all interpretations allowed.
					baseval = 0.5;
				}
			}

			if (!baseval)
			{
				continue;
			}

			// There should be some other flag for the actual search depth
			if (genwidth == HAPLOTYPING)
			{
				if (zeropropagate && extparams.gstr)
				{
					*(extparams.gstr) *= 2;
					*(extparams.gstr) += strain - 1;
				}
			}
			else
			{
				// Track to the parent generation, creating evaluation objects first
				// These do a lookup in a special hash for combinations known to be 0, to avoid unnecessary calls
				// Both are checked before either branch of the pedigree is actually traced.
				recursetrackpossible<update, zeropropagate> subtrack1 =
					recursetrackpossible<update, zeropropagate>(this, tb, markerval, marker,
					upflag,
					upflag2,
					upshift,
					genwidth,
					f2n,
					firstpar,
					0,
					extparams);

				if (subtrack1.prelok && (!zeropropagate || (genwidth == 1 << (NUMGEN - 1))) )
					baseval *= recursetrackpossible<update, zeropropagate>(this, tb, themarker[!realf2n], marker,
					upflag,
					upflag2,
					upshift,
					genwidth,
					f2n,
					!firstpar,
					1,
					extparams);
				if (!baseval) continue;

				baseval *= subtrack1;
			}

			if (!baseval) continue;

			ok += baseval;

			if (update && !allthesame)
			{
				(*tb.haplos)[n][f2n] += extparams.updateval;
			}
		}

		return ok;
	}


	// calltrackpossible is a slight wrapper that hides at least some of the interanl parameters needed for the recursion from outside callers
	template<bool update, bool zeropropagate> double calltrackpossible(const threadblock& tb, const int* markervals, const unsigned int marker,
		const int genotype, const unsigned int offset, const int flag2, const double updateval = 0.0)
	{		
		return trackpossible<update, zeropropagate>(tb, 0, marker, genotype * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(updateval, 0));
	}

	// "Fix" parents, i.e. infer correct marker values for any zero values existing.
	// Depending on some specifics, this can either assume that all existing marker values will be reflected in some offspring individuals,
	// or indeed only infer specifically those things that can be said from certain, i.e. if a zero value exists, and the other branch of
	// the pedigree doesn't allow a marker value found in the offspring, then it can be assumed to have been transmitted from the individual
	// with the zero value. 
	void fixparents(unsigned int marker, bool latephase)
	{
		pair<int, int>& themarker = markerdata[marker];
		individ* parp[3] = {pars[0], pars[1], this};
		int az = 0;

#pragma omp critical (parmarkerval)
		for (int i = 0; i < 3; i++)
		{
			if (parp[i])
			{
				az += parp[i]->markerdata[marker].first == 0;
				az += parp[i]->markerdata[marker].second == 0;
				while (parp[i]->markervals.size() <= marker)
				{
					parp[i]->markervals.resize(markerdata.size());
				}
			}
		}

		//if (!az) return;		

		double okvals[2] = {0};
		// We only accept an interpretation when it is by exclusion the only possible one. As soon as one intepretation has gained acceptance,
		// no need to retest it.
		for (shiftflagmode = 0; shiftflagmode < 1; shiftflagmode+=1)
		{
		  			threadblock tb;
			for (int i = 0; i < NUMTYPES; i++)
			{
				for (int flag2 = 0; flag2 < NUMPATHS; flag2++)
				{
					if (okvals[flag2 & 1]) continue;

					double ok = calltrackpossible<false, false>(tb, 0, marker, i, 0, flag2);
					if (!ok) continue;
					okvals[flag2 & 1] += ok;
					if (okvals[0] && okvals[1]) break;
				}
				if (okvals[0] && okvals[1]) break;
			}
		}

		if (!okvals[0] && !okvals[1])
		{
			printf("Clearing %d:%d\n", this->n, marker);
			markerdata[marker] = make_pair(0, 0);
		}

		if ((((bool) okvals[0]) ^ ((bool) okvals[1])) || latephase)
		{
			for (int flag2 = 0; flag2 < 2; flag2++)
			{
				if (!okvals[flag2]) continue;

				for (int k = 0; k < 2; k++)
				{
					if (pars[k])
					{
					  int u = ((k ^ flag2 /*^ *tb.shiftflagmode*/) & 1);
						if ((&themarker.first)[u])
				  {
#pragma omp critical (parmarkerval)
					  pars[k]->markervals[marker].insert((&themarker.first)[u]);
				  }
					}
				}
			}
		}		
	}

	// Verifies that an update should take place and performs it. Just a wrapper to calltrackpossible.
	void updatehaplo(const threadblock& tb, const unsigned int marker, const unsigned int i, const int flag2, const double updateval)
	{
		pair<int, int>& themarker = markerdata[marker];		

		double ok = calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);

		if (ok)
		{
			calltrackpossible<true, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);
		}
		else
		{

		}
	}

	// Adjust the probability, i.e. filter all probability values based on the haplotype weights and overall admissibility for the different
	// states.
	void adjustprobs(const threadblock& tb, double* probs, const unsigned int marker, double& factor, const bool oldruleout, int flag99)
	{
		double sum = 0;
		double probs2[NUMTYPES];

		const pair<int, int>& themarker = markerdata[marker];
		const bool ruleout = true; // TODO

		for (int q = 0; q <= (int) !ruleout; q++)
		{
			int f2start = 0;
			int f2end = NUMPATHS;

			// Negative values other than -1 just cause trouble

			// Can we really trust the logic to expand correctly even for zero values?
			//if (flag99 >= 0 || !somezero[marker])
			{
				f2start = flag99;
				f2end = flag99 + 1;
			}

			for (unsigned int i = 0; i < NUMTYPES; i++) // genotype
			{
				probs2[i] = probs[i];

				// We will multiply this already small number with an even smaller number... let's assume it's zero and be done with it.
				if (probs[i] < 1e-30)
				{
					probs[i] = 0;
					continue;
				}

				double realok = 0;

				for (int flag2 = f2start; flag2 < f2end; flag2++)
				{
					realok += calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2);
				}					

				// TODO UGLY CAPPING
				if (HAPLOTYPING)
				{
					probs[i] *= (double) realok;
				}
				else
				{
					probs[i] *= (bool) realok;
				}
				sum += probs[i];
			}

			if (sum == 0 && !ruleout)
			{
				for (int i = 0; i < NUMTYPES; i++)
				{
					probs[i] = probs2[i];
				}

				// The code sees: This doesn't make sense, ignore this marker!
				// NOTE: If the parent-grandparent genotypes are completely incompatible, the problem
				// is NOT eliminated.

				// Current caching scheme does not execute full elimination. Efficiency gains possible by not allowing
				// for genotyping error in every single marker. A small epsilon should otherwise be introduced whenever
				// any probs[i] is zero.
				markerdata[marker] = make_pair(0, 0);
				fprintf(stderr, "Error in %x, marker %d, impossible path\n", this, marker);
				sum = 1;
			}
			else
				break;
		}

		// Normalize, update auxiliary exponent
		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] /= sum;
		}
		factor += log(sum);
	}

	// Append a "multi-step" transition. If the cached values (essentially a N * N transition matrix for the steps from startmark to
	// endmark) are missing, calculate them first.
	double fillortake(const threadblock& tb, const int index, const unsigned int startmark, const unsigned int endmark, double* probs)
	{
		if ((tb.done[*(tb.shiftflagmode)])[index] != (*tb.generation))
		{

			for (unsigned int i = 0; i < NUMTYPES; i++)
			{
				double probs2[NUMTYPES] = {0};
				probs2[i] = 1;

				// Note ruleout here, could be set to false if we "prime" with an evenly distributed probs first
				(tb.factors[*tb.shiftflagmode])[index][i] = quickanalyze<false, noneturner>(tb, none, startmark,
					endmark,
					-1,
					0,
					-1,
					true,
					probs2);

				double sum = 0;

#pragma ivdep:back
				for (int j = 0; j < NUMTYPES; j++)
				{
					if (!_finite(probs2[j])) probs2[j] = 0.0;
					sum += probs2[j];
				}

				sum *= NUMTYPES;

				(tb.factors[*tb.shiftflagmode])[index][i] += log(sum);
				float* probdest = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];

#pragma ivdep
				for (int j = 0; j < NUMTYPES; j++)
				{
					probdest[j] = probs2[j] / sum;
				}
			}
			(tb.done[*tb.shiftflagmode])[index] = (*tb.generation);
		}

		float factor = MINFACTOR;
		for (int i = 0; i < NUMTYPES; i++)
		{
			factor = max(factor, (tb.factors[*tb.shiftflagmode])[index][i]);
		}

		double probs2[NUMTYPES] = {0};

		for (int i = 0; i < NUMTYPES; i++)
		{
			float step = (tb.factors[*tb.shiftflagmode])[index][i] - factor;
			if (probs[i] == 0.0 || step < -20.0f) continue;
			double basef = exp((double) step) * probs[i];
			//			if (basef == 0.0 || !_finite(basef)) continue;
			const float* probsource = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];
#pragma ivdep
			for (int j = 0; j < NUMTYPES; j++)
			{
				const double term = basef * probsource[j];
				probs2[j] += term;
			}
		}

		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] = probs2[i];
		}


		return factor;
	}

	// Is it OK to cache over this range, i.e. is no fixed position found here. For now, lockpos is an integer, but it can be changed to a vector
	// for actually fixing multiple loci.
	const bool okstep(const int startmark, const int endmark, const int lockpos) const
	{
		bool found = (lockpos <= -1000 - startmark && lockpos > -1000 - endmark) || (lockpos >= markerposes[startmark] && lockpos <= markerposes[endmark]);

		return !found;
	}

	// Analyze for a specific range, including a possible fixed specific state at some position (determined by lockpos and genotype)
	template<bool inclusive, class T> double quickanalyze(const threadblock& tb, const T& turner, unsigned int startmark,
		const unsigned int endmark, const int lockpos, const int genotype, const int flag2, bool ruleout, double* probs,
		float minfactor = MINFACTOR)
	{
		unsigned int stepsize;
		double factor = 0;
		bool allowfull = inclusive;
		bool frommem = false;
		if (inclusive && tb.quickgen[*tb.shiftflagmode] == *tb.generation && tb.lockpos[*tb.shiftflagmode] == lockpos)
		{
			allowfull = true;
			frommem = true;
			factor = tb.quickfactor[*tb.shiftflagmode];
			if (factor <= minfactor)
			{
				return MINFACTOR;
			}

			startmark = tb.quickmark[*tb.shiftflagmode];
			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] = tb.quickmem[*tb.shiftflagmode][i];
			}
		}

		// Loop until we have reached the end, with blocks of varying sizes.
		while (startmark < endmark)
		{
			for (stepsize = 1; stepsize < (endmark - startmark + allowfull) &&
				okstep(startmark, startmark + stepsize, lockpos) &&
				!(startmark & (stepsize - 1)); stepsize *= 2);

			// A single step, either due to the fixated genotypee being within this range, or simply because we've gone all the way down
			// the tree.
			if (stepsize <= 2)
			{
				stepsize = 1;

				if (!frommem && !okstep(startmark, startmark + 1, lockpos))
				{
					// If we have a fixated genotype at marker x in one call, it is very likely that the next call will also
					// be a fixated genotype at marker x. Possibly another one, but still fixated. The fixated genotype does not
					// change the probability values leading up to this position, so they are cached.
					tb.quickgen[*tb.shiftflagmode] = *tb.generation;
					tb.quickmark[*tb.shiftflagmode] = startmark;
					tb.lockpos[*tb.shiftflagmode] = lockpos;
					tb.quickfactor[*tb.shiftflagmode] = factor;

					for (int i = 0; i < NUMTYPES; i++)
					{
						tb.quickmem[*tb.shiftflagmode][i] = probs[i];
						tb.quickendfactor[*tb.shiftflagmode][i] = 1.0;
					}
					frommem = true;
				}
				bool willquickend = (frommem && lockpos <= -1000 && genotype >= 0);

				if (willquickend)
				{
					if (factor + tb.quickendfactor[*tb.shiftflagmode][genotype] <= minfactor)
					{
						return MINFACTOR;
					}
					// If we are doing a quick end
					factor += realanalyze<4, T>(tb, turner, startmark, startmark + stepsize, lockpos, genotype, flag2, ruleout, probs);
				}
				else
				{
					factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, lockpos, genotype, flag2, ruleout, probs);
				}

				// This will work.
				if (!_finite(factor) || factor <= minfactor)
				{
					//			    			    if (startmark) printf("%d;%d;%d : %lf\n", n, shiftflagmode, startmark, factor);
					return MINFACTOR;
				}
				if (willquickend)
				{
					double wasval = probs[genotype];
					if (tb.quickendfactor[*tb.shiftflagmode][genotype] > 0.0)
					{
						// Precalc the full set of transitions from this very marker, for all states, all the way to the end
						// That way, any new analysis starting from this marker can be solved with no recomputation at all.
						// Note that only the latest marker used is stored here, but as we generally loop over all states and possibly
						// multiple flag values, it still helps.
						for (int i = 0; i < NUMTYPES; i++)
						{
							tb.quickendprobs[*tb.shiftflagmode][genotype][i] = 0;				     
						}
						tb.quickendprobs[*tb.shiftflagmode][genotype][genotype] = 1.0;				 
						double superfactor = realanalyze<0 | 2, T>(tb, turner, startmark, startmark + stepsize, -1, -1, -1, ruleout, &tb.quickendprobs[*tb.shiftflagmode][genotype][0]);

						superfactor += quickanalyze<true, T>(tb, turner, startmark + stepsize, endmark, -1, -1, flag2, ruleout, &tb.quickendprobs[*tb.shiftflagmode][genotype][0],
							// all uses of this precalced data will have the
							// quickfactor component in common, so taking that
							// into account for the limit is not problematic
							// any strong filtering at this specific locus can be
							// handled by the return MINFACTOR line above
							minfactor - tb.quickfactor[*tb.shiftflagmode]);

						tb.quickendfactor[*tb.shiftflagmode][genotype] = superfactor;
					}

					factor += tb.quickendfactor[*tb.shiftflagmode][genotype];
					factor += log(wasval);
					if (factor <= minfactor) return MINFACTOR;

					for (int i = 0; i < NUMTYPES; i++)
					{
						probs[i] = tb.quickendprobs[*tb.shiftflagmode][genotype][i];
					}

					return factor;
				}
			}
			else
			{
				// Use a stored multi-step transition.
				stepsize /= 2;

				int base = 0;
				for (int q = stepsize / 2; q > 1; q /= 2)
				{
					base += markerposes.size() / q;
				}

				base += startmark / stepsize;

				factor += fillortake(tb, base, startmark, startmark + stepsize, probs);
			}
			startmark += stepsize;
			allowfull |= true;

			if (!_finite(factor) || factor <= minfactor)
			{
				if (!frommem && !okstep(startmark, endmark, lockpos))
				{
					tb.quickgen[*tb.shiftflagmode] = *tb.generation;
					tb.lockpos[*tb.shiftflagmode] = lockpos;
					tb.quickfactor[*tb.shiftflagmode] = MINFACTOR;
				}

				return MINFACTOR;
			}
		}		

		return factor;
	}

	// A wrapper to quickanalyze, preparing the start probability vector.
	template<class T> double doanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const int lockpos,
		const int genotype, const int flag2, bool ruleout = false, double* realprobs = 0, float minfactor = MINFACTOR)
	{
		double fakeprobs[NUMTYPES];
		double* probs;

		if (realprobs)
		{
			probs = realprobs;
		}
		else
		{
			probs = fakeprobs;
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN;
			}
		}

		double factor = quickanalyze<true, T>(tb, turner, startmark, endmark, lockpos, genotype, flag2, ruleout, probs, minfactor);
		bool small = !_finite(factor) || minfactor >= factor;

		if (!small) adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO, the very last marker can be somewhat distorted

		return factor;
	}

	// This is the actual analyzing code. It works with no caches, and can function independently, but is generally only used to patch in those
	// elements not found in caches by quickanalyze and fillortake.
	//
	// Both transition and emission (through adjustprobs) probabilities are handled here.
	//
	// first bit in updateend signals whether the interval is end-inclusive at endmark
	// the second bit in updateend will be 0 if the interval is end-inclusive at startmark, and 1 IF NOT
	// the third bit will cause the code to quit early, after processing the genotype and turner condition
	template<int updateend, class T> double realanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const int lockpos, const int genotype,
		const int flag2, const bool ruleout = false, double* realprobs = 0)
	{
		double fakeprobs[NUMTYPES];
		double* probs;

		if (realprobs)
		{
			probs = realprobs;
		}
		else
		{
			probs = fakeprobs;
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN;
			}
		}

		// The *logged* normalization factor
		double factor = 0;

		// Walk over all markers.
		for (int j = startmark + 1; j <= endmark; j++)
		{
			double startpos = markerposes[j - 1];
			double endpos = markerposes[j];

			bool tofind = (lockpos >= startpos && lockpos <= endpos);
			if (tofind)
			{
				endpos = lockpos;
			}

			int f2use = -1;

			if (lockpos == -1000 - j + 1)
			{
				tofind = true;

				endpos = startpos;

				f2use = flag2;
			}

			// If we are at the very first position, and the specific flag was set, include the emission probabilities for the previous
			// marker. Used to maximize the caching.
			if (!((updateend & 2) && (j == startmark + 1))) adjustprobs(tb, probs, j - 1, factor, ruleout, f2use);

			// For a specific intra-marker region, we have two cases: the case of a fixated position between the two markers, and the simple case
			// of no fixated position.
			for (int iter = 0; iter <= (int) tofind; iter++)
			{
				// If iter is 1, we have currently handled the transition all the way to the fixated position. Now filter to keep only
				// a single state value positive.
				if (iter)
				{
					turner(probs);
					if (genotype >= 0)
					{
						for (int i = 0; i < NUMTYPES; i++)
						{
							probs[i] = probs[i] * (i == genotype);
						}
					}

					// Were we asked to stop at this very state, in the middle of things?
					if (updateend & 4) return factor;
				}


				double dist = (endpos - startpos);

				// Compute transitions, assuming there is any distance to transfer over.
				if (dist > 0)
				{
					double probs2[NUMTYPES] = {0};
					double recprob[2];

#pragma ivdep
					// Compute recombination probabilities for this specific distance, for the two sexes.
					// (as the sex-dependent marker distance might not be a simple transformation, the actrec
					// data comes into play).
					for (int k = 0; k < 2; k++)
					{
						recprob[k] = 0.5 * (1.0 - exp(actrec[k][j] * (dist)));
						//					if (iter == tofind) recprob[k] = max(recprob[k], 1e-5);
					}

					// Precompute values for presence (/lack of) a crossover event for either sex.
					double other[2][2];
#pragma ivdep
					for (int m = 0; m < 2; m++)
					{
						for (int k = 0; k < 2; k++)
						{
							double prob = recprob[k];
							if (m) prob = 1.0 - prob;

							other[m][k] = prob;
						}
					}

					double recombprec[NUMTYPES];
#pragma ivdep
					for (int index = 0; index < NUMTYPES; index++)
					{
						recombprec[index] = 1;
					}


					// Compute probabilities for arbitrary xor values of current and future state
					for (int t = 0; t < TYPEBITS; t++)
					{
						int sex = TYPESEXES[t];

#pragma ivdep
						for (int index = 0; index < NUMTYPES; index++)
						{
							int val = !((index >> t) & 1);
							recombprec[index] *= other[val][sex];
						}

					}

					// Use those xor values
					// For the 4-state model, this is an inefficient way to go about it, but it is quite a bit more efficient for
					// the 64-state model (or beyond).
					for (int from = 0; from < NUMTYPES; from++)
					{
						if (probs[from] < MINFACTOR) continue;
						for (int to = 0; to < NUMTYPES; to++)
						{
							probs2[to] += probs[from] * recombprec[from ^ to];
						}
					}

					for (int c = 0; c < NUMTYPES; c++)
					{
						probs[c] = probs2[c];
					}
				}
				//				else
       				{
				  //				if (iter == tofind)
					{
						for (int c = 0; c < NUMTYPES; c++)
						{
							if (probs[c] < 1e-7) probs[c] = 1e-7;
						}
					}
				}

				startpos = endpos;
				endpos = markerposes[j];
			}			
		}

		if (updateend & 1)
		{
			adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO
		}

		return factor;
	}
};

// Oh, how we waste memory, in a pseudo-O(1) manner
individ* individer[2000000];
// dous contains those individuals that really should be analyzed
vector<individ*> dous;

// retrieve the individual with a specific number
individ* const getind(int n)
{
	if (n <= 0) return 0;

	if (!individer[n])
	{
		individer[n] = new individ();
		individer[n]->n = n;
	}

	return individer[n];
}

// read qtlmas sample data, code used for testing haplotypes, quite a bit of hardcoding present in this function
void readqtlmas()
{
	FILE* inpheno = fopen("phenotype.txt", "r");
	FILE* ingeno = fopen("genotype_cor.txt", "r");
	markertranslation.resize(6000);
	for (int i = 0; i < 6; i++)
	{
		chromstarts.push_back(i * 1000);
		for (int j = 0; j < 1000; j++)
		{
			markertranslation[i * 1000 + j] = i * 1000 + j;
			markerposes.push_back(j * 0.1);
			for (int t = 0; t < 2; t++)
			{
				actrec[t].push_back(baserec[t]);
			}
		}
	}
	chromstarts.push_back(6000);

	char tlf[16384];
	fgets(tlf, 16384, inpheno);

	for (int indn = 1; indn < 7000; indn++)
	{
		int fa = 0;
		int mo = 0;
		individ* ind = getind(indn);
		ind->pars[0] = getind(fa);
		ind->pars[1] = getind(mo);

		ind->gen = 5;

		ind->sex = 0;
		ind->strain = 1;

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());
		ind->markerdata.resize(markerposes.size());

		for (int i = 0; i < 6000; i++)
		{
			ind->haploweight[i] = 0.5;
		}
	}

	while (fgets(tlf, 16384, inpheno))
	{
		int indn, fa, mo, sex, gen;
		if (sscanf(tlf, "%d %d %d %d %d", &indn, &fa, &mo, &sex, &gen) < 4) break;


		individ* ind = getind(indn);

		ind->pars[0] = getind(fa);
		ind->pars[1] = getind(mo);


		if (fa && mo && (!ind->pars[0]->haplobase.size() ||
			!ind->pars[1]->haplobase.size()))
		{
			printf("PROBLEM %d %d %d\n", indn, fa, mo);
		}
		ind->gen = gen;

		ind->sex = sex - 1;
		ind->strain = (indn % 2) + 1;

		if (gen >= 1) dous.push_back(ind);

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());
		ind->markerdata.resize(markerposes.size());
	}

	int indn;
	while (fscanf(ingeno, "%d", &indn) == 1 && indn)
	{
		individ* ind = getind(indn);
		if (!ind->markerdata.size())
		{
			printf("Stopping at %d\n", indn);
			break;
		}

		for (int i = 0; i < 6000; i++)
		{
			int a, b;
			fscanf(ingeno, "%d %d", &a, &b);

			ind->markerdata[i] = make_pair(a, b);

		}
	}
}


// read marker info in a format similar to that accepted by ccoeff
void readmarkerinfo(FILE* in)
{
	int n, m;
	fscanf(in, "%d %d", &n, &m);
	vector<int> count;
	markertranslation.resize(m);

	for (int i = 0, j = 0; i < n; i++)
	{
		fscanf(in, "%d", &m);
		count.push_back(m);

		printf("Number on ch %d: %d \n", i + 1, m);

		for (int k = 0; k < count.end()[-1]; k++)
		{
			fscanf(in, "%d", &m);
			markertranslation[m - 1] = ++j;
		}
	}

	int pos = 0;
	for (int i = 0; i < n; i++)
	{
		chromstarts.push_back(pos);		

		vector<double> part[2];
		for (int t = 0; t < sexc; t++)
		{
			fscanf(in, "%d", &m);

			double j = 0;
			for (int k = 0; k < count[i]; k++)
			{
				double v;
				fscanf(in, "%lf", &v);
				j += v;

				for (int p = 0; p < 2 / sexc; p++)
				{
					part[t + p].push_back(j / discstep);
				}
			}
		}

		for (int k = 0; k < count[i]; k++)
		{
			double sum = 0;
			for (int t = 0; t < 2; t++)
			{
				sum += part[t][k];
			}
			sum /= 2.0;

			markerposes.push_back(sum);
			printf("%lf\n", sum);

			for (int t = 0; t < 2; t++)
			{
				double avgdist = 0;
				if (k)
				{
					avgdist = markerposes[pos] - markerposes[pos - 1];
				}

				if (avgdist)
				{					
					actrec[t].push_back(baserec[t] * (part[t][k] - part[t][k - 1]) / avgdist);
				}
				else
				{
					// This should give a loud error.
					actrec[t].push_back(-1.0);
				}
			}
			pos++;
		}
	}
	chromstarts.push_back(pos);

	printf("%d chromosomes, %d markers, really %d/%d\n", n, m, markerposes.size(), chromstarts.size() - 1);
}

// read a pedigree in a format similar to the one accepted by ccoeff
void readped(FILE* in)
{
	int famsize;

	// The format is a set of full sibships, each started by four founders (possibly some identical), two parents
	while (fscanf(in, "%d", &famsize) == 1)
	{
		for (int i = 0; i < famsize + 6; i++)
		{
			int indn;
			int fa, mo;
			int sex = 0;

			// strain here could just as well have been called line
			// OTOH, in a line-based file format, animal lines can be confusing
			int strain = -1;
			char tlf[255];
			while (fgets(tlf, 255, in))
			{
				if (sscanf(tlf, "%d %d %d %d %d", &indn, &fa, &mo, &sex, &strain) >= 3) break;
			}

			individ* ind = getind(indn);
			ind->pars[1] = getind(fa);
			ind->pars[0] = getind(mo);

			if (ind->pars[0] && ind->pars[1] && ind->pars[0]->sex == 1 && ind->pars[1]->sex == 0) {
			  swap(ind->pars[0], ind->pars[1]);
			  printf("%d\n", indn);
			}

			ind->sex = sex - 1;
			ind->strain = strain;
			// The generation is defined by the index, the first few are always the first two generations
			ind->gen = 0;
			ind->gen += i >= 2;
			ind->gen += i >= 6;

			// Only analyze the F2s
			if (i >= 6) dous.push_back(ind);
		}
	}
	printf("Pedigree containing %d F2 individuals\n", dous.size());
}

// Read marker data and dimension several data fields dependent on marker data
// These are expected to be in a format similar to the one accepted by ccoeff
void readmarkerdata(FILE* in)
{
	int indn;
	int n = 0;
	while (fscanf(in, "%d", &indn) == 1)
	{
		individ* ind = getind(indn);
		n++;
		ind->markerdata.resize(markerposes.size());
		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		for (unsigned int i = 0; i < markerposes.size(); i++)
		{
			ind->haploweight[i] = 0.5;
		}

		for (unsigned int i = 0; i < markertranslation.size(); i++)
		{
			int a, b;
			if (fscanf(in, "%d %d", &a, &b) != 2)
			{
				fprintf(stderr, "Marker mismatch: %d\n", indn);
			}
			if (markertranslation[i])
			{
				ind->markerdata[markertranslation[i] - 1] = make_pair(a,b);
			}
		}
	}
	printf("Marker data parsed for %d individuals\n", n);
}

// Some operations performed when marker data has been read, independent of format.
void postmarkerdata()
{
	int any, anyrem;
	bool latephase = false;

	// If inference is active, add "new" marker data until all data has been found.
	if (CORRECTIONINFERENCE) do
	{
#pragma omp parallel for schedule(dynamic,32)
		// all haploweights MUST be non-zero at this point, as we do not explore all shiftflagmode values
		for (int i = 1; i < 2000000; i++)
		{
			individ* ind = getind(i);
			if (ind->markerdata.size())
			{
				generation++;
				
				for (int g = 0; g < ind->markerdata.size(); g++)
				{
					ind->fixparents(g, latephase);
				}
			}
		}

		any = 0;
		anyrem = 0;
		for (int i = 1; i < 2000000; i++)
		{
			individ* ind = getind(i);

			for (int g = 0; g < (int) ind->markervals.size(); g++)
			{
				const int startsize = ind->markervals[g].size();
				const int known =
					((bool) ind->markerdata[g].first) +
					((bool) ind->markerdata[g].second);

				const int oldany = any;

				if (latephase && known == 2 && !startsize && ind->gen < 2)
				{
					ind->markerdata[g].second = 0;
					ind->markerdata[g].first = 0;
					any++;
					anyrem++;
				}
				else
					if (known == 2) continue;

				if (ind->markerdata[g].first) ind->markervals[g].insert(ind->markerdata[g].first);
				if (ind->markerdata[g].second) ind->markervals[g].insert(ind->markerdata[g].second);

				if (ind->markervals[g].size() >= 3)
				{
					fprintf(stderr, "Error, too many matches: %d\t%d\n", i, g);
				}
				if (!latephase && ind->markervals[g].size() == 2)
				{
					ind->markerdata[g] = make_pair(*ind->markervals[g].begin(), *(++ind->markervals[g].begin()));
					any++;
				}
				if (latephase && ind->markervals[g].size() == 1 && startsize == 1 && known == 1 && !ind->pars[0] && !ind->pars[1])
				{
					if (!ind->markerdata[g].first || !ind->markerdata[g].second) any++;
					ind->markerdata[g] = make_pair(*ind->markervals[g].begin(), *ind->markervals[g].begin());					
				} // DANGEROUS ASSUMPTIONS
				else if (!latephase && ind->markervals[g].size() == 1 && known == 0)
				{
					any++;
					ind->markerdata[g] = make_pair(*ind->markervals[g].begin(), 0);
				}

				if (any != oldany) printf("Correction at %d, marker %d (%d;%d)\n", i, g,
					ind->markerdata[g].first, ind->markerdata[g].second);

			}

			for (int g = 0; g < (int) ind->markervals.size(); g++)
			{
			  if (ind->markerdata[g].first == sexmarkerval) {
			    ind->markerdata[g] = make_pair(ind->markerdata[g].second, ind->markerdata[g].first);
			  }
				ind->markervals[g].clear();
			}
		}
		fprintf(stderr, "Number of corrected genotypes: %d\n", any);
		if (latephase)
		{
			latephase = false;
		}
		else
		{
			if (!any && !latephase)
			{
				any++;
				latephase = true;
			}
		}
	}
	while (any > anyrem);

	for (int i = 1; i < 2000000; i++)
	{
		individ* ind = getind(i);
		ind->markervals.clear();

		// Lock the first position
		if (HAPLOTYPING && ind && ind->haploweight.size())
			// These days we do the locking in all generations
		{
			// Lock the haplotype (in an arbitrary manner) for the first marker in each linkage group
			// Propagation would be more efficient if we locked the mid-position (which should really be determined in cM)
			// Doing so would be much more opaque, though...
			for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
			{
				unsigned int j;
				for (j = chromstarts[i]; j != chromstarts[i + 1] && ind->markerdata[j].first == ind->markerdata[j].second; j++);

				if (j != chromstarts[i + 1]) ind->haploweight[j] = 0;
			}
		}
	}
}

// The actual walking over all chromosomes for all individuals in "dous"
// If "full" is set to false, we assume that haplotype inference should be done, over marker positions.
// A full scan is thus not the iteration that takes the most time, but the scan that goes over the full genome grid, not only
// marker positions.
template<bool full> void doit(FILE* out)
{
	const bool doprint = full || true;

	int count = 0;
	for (unsigned int j = 0; j < dous.size(); j++)
	{
		if (dous[j]->markerdata.size())
		{
			count++;
		}
		else
		{
			dous.erase(dous.begin() + j);
		}
	}


	if (doprint)
	{
		fprintf(out, "%d %d\n", count, chromstarts.size() - 1);
	}

	for (int i = 0; i < 2000000; i++)
	{
		individ* ind = getind(i);
		if (!ind) continue;

		for (unsigned int j = 0; j < ind->haplocount.size(); j++)
		{
			ind->haplocount[j] = 0.0f;
			ind->haplobase[j] = 0.0f;
		}
	}


	for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
	{
		if (doprint)
		{
			fprintf(out, "%d %d\n", i + 1, (int) markerposes[chromstarts[i + 1] - 1]);
		}
		printf("Chromosome %d\n", i + 1);

		// The output for all individuals in a specific iteration is stored, as we have parallelized the logic and 
		// want the output to be done in order.
		vector<array<char, 150000> > outqueue;
		vector<int> oqp; // output queue position

		oqp.resize(dous.size());
		outqueue.resize(dous.size());


#pragma omp parallel for schedule(dynamic,1)
		for (int j = 0; j < (int) dous.size(); j++)
		{
			generation++;
			threadblock tb;

			// Some heaps are not properly synchronized. Putting a critical section here makes the operations not safe,
			// but *safer*.
#pragma omp critical(uglynewhack)
			for (int t = 0; t < NUMSHIFTS; t++)
			{
				factors[t].resize(markerposes.size());
				cacheprobs[t].resize(markerposes.size());
				done[t].resize(markerposes.size());
			}

			if (dous[j]->markerdata.size())
			{				
				int qstart = -1000 - chromstarts[i];
				int qend = -1000 - chromstarts[i + 1];
				int qd = -1;
				int f2s = 0;
				int f2end = NUMPATHS;

				if (!HAPLOTYPING)
				{
					f2s = -1;
					f2end = 0;
				}

				int shifts = 0;
				int shiftend = NUMSHIFTS;

				reltree.clear();
				reltree.push_back(dous[j]);
				int flag2ignore = 0;

				// Special optimization hardcoded for this population structure, eagerly skipping flags that do not
				// correspond to any inheritance, i.e. if not the full pedigree of 6 individuals back is present.
				//				if (HAPLOTYPING && NUMGEN == 3)
				{
					flag2ignore = 1;
					for (int lev1 = 0; lev1 < 2; lev1++)
					{
						individ* lev1i = dous[j]->pars[lev1];
						if (!lev1i) continue;
						int flag2base = 1 << (1 + lev1 * NUMFLAG2GEN);
						flag2ignore |= flag2base;

						reltree.push_back(lev1i);
						for (int lev2 = 0; lev2 < 2; lev2++)
						{
							individ* lev2i = lev1i->pars[lev2];
							if (!lev2i) continue;

							flag2ignore |= (flag2base << (lev2 + 1));
							reltree.push_back(lev2i);
						}
					}

					flag2ignore ^= (NUMPATHS - 1);
				}

				sort(reltree.begin(), reltree.end());
				reltree.resize(unique(reltree.begin(), reltree.end()) - reltree.begin());

			        bool skipsome = false;
				for (int u = 0; u < reltree.size() && !skipsome; u++)
				  {
				    for (int p = 0; p < 2; p++)
				      {
					if ((&(reltree[u]->markerdata[chromstarts[i]].first))[p] == sexmarkerval)
					  skipsome = true;
				      }
				  }

				if (!(HAPLOTYPING && NUMGEN == 3))
				  {
				    reltree.resize(0);
				    flag2ignore = 0;
				  }

				if (full)
				{
					qstart = (int) markerposes[chromstarts[i]];
					qend = (int) markerposes[chromstarts[i + 1] - 1] + 1;
					qd = 1;
					f2s = -1;
					f2end = 0;
					/*shifts = 0;
					shiftend = 1;*/
				}

				if (dous[j]->gen < 2) shiftend = 2;

				double factor = -1e15;
				double factors[NUMSHIFTS];
				for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
				{
					factors[shiftflagmode] = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, -1, 0, -1, false, 0, -20 + factor);
					factor = max(factor, factors[shiftflagmode]);
				}

				// Normalize!
				double realfactor = 0;
				for (int s = shifts; s < shiftend; s++)
				{
					realfactor += exp(factors[s] - factor);
				}


				printf("%d,%03d: %lf\t", dous[j]->n, flag2ignore, factor);
				factor += log(realfactor);
				printf("%lf\n", factor);
				//fflush(stdout);
				// Flushing can be useful for debugging, but not for performance!
				// This output can get ugly due to race conditions. One shouldn't rely on it.

				if (_isnan(factor)) continue;

				char lineout[255];

				// States are mapped onto values describing the line/strain origin, in the sense of 00, 01, 10 or 11
				int maptogeno[NUMTYPES];
				shiftflagmode = 0;
				for (int g = 0; g < NUMTYPES; g++)
				{
					int sum = 0;
					dous[j]->trackpossible<false, true>(tb, 0, 0, g * 2, 0, 0, trackpossibleparams(0, &sum));

					// switch due to earlier inner switching
					if (sum == 2 || sum == 1) sum = 3 - sum;

					maptogeno[g] = sum;
				}

				// Walk over all chromosome positions, whether it be markers (negative q values <= -1000) or grid positions
				for (int q = qstart; q != qend; q+=qd)
				{
					double probs[4] = {0};

					for (int g = 0; g < NUMTYPES; g++)
					{						
						for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
						{
							if (factor - factors[shiftflagmode] > 20) continue;
							for (int flag2 = f2s; flag2 < f2end; flag2++)
							{
								if (flag2 & (flag2ignore)) continue;

								int firstpar = 0;
								double val;

								if (q <= -1000)
								{
									// Do a lookup in the impossible structure. If we are at a marker position, those branches that are impossible will have an exact
									// score of 0.
									// If we are NOT at a marker, this information can still be used in trackpossible, but not in here, as there is some recombination
									// added between the marker and fixated position we're checking.
									for (int a = 0; a < 2; a++)
									{
										if (a) firstpar = !firstpar;
										const int genwidth = (1 << (NUMGEN - 1));

										int f2n = ((flag2 ^ shiftflagmode) & 1);

										int upflagr = upflagit(g, firstpar, genwidth);
										int upflag2r = upflagit(flag2 >> 1, firstpar, genwidth);
										int upshiftr = upflagit(shiftflagmode >> 1, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));
										int marker = -(q + 1000);
										int impossibleval = generation * markerposes.size() + marker;

										if (impossible[shiftflagmode & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3] == impossibleval)
										{
											goto continueloop;
										}
									}
								}

								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, q, g,
									flag2, true, 0, -20.0 + factor) - factor;

								if (_finite(val) && val > -20.0)
								{
									val = exp(val);
									int mapval = maptogeno[g];
									for (int p = 0; p < 2; p++)
									  {
									    if ((&(dous[j]->markerdata[chromstarts[i]].first))[p] == sexmarkerval)
									      {
										mapval &= 3 - (1 << ((shiftflagmode ^ p /*^ 1*/) & 1));
									      }
									  }
									probs[mapval] += val;
									//									printf("%d %d %d %d %lf\n", j, q, g, flag2, val);
									if (!full) dous[j]->updatehaplo(tb, -q - 1000, g, flag2, val);
								}
continueloop:;
							}
						}
					}

					// Consider doing haplotype reversal from a specific position and all the way down.
					if (!early && !full && dous[j]->gen >= 1)
					{
						double rawvals[NUMTYPES * 2][NUMSHIFTS];
						double sumnegval[TYPEBITS + 1] = {0};
						for (int g = 0; g < NUMTYPES * 2; g++)
						{
							for (int s = 0; s < NUMSHIFTS; s++)
							{
								rawvals[g][s] = 0;
							}
						}

						for (int g = 0; g < NUMTYPES * 2; g++)
						{		
							if (g & (flag2ignore >> 1)) continue;

							int c = 0;
							for (int p = 0; p < TYPEBITS + 1; p++)
							{
								if (g & (1 << p)) c++;
							}

							if (c > 1) continue;

							aroundturner turn(g);
							for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
							{
								// If we are above this limit, we are shifting shift mode
								// within the range and cannot use this heuristic of the
								// aggregated probability to know anything... anything at all!
								// (In other cases, our simulated crossover event amounts to 
								// far below a exp(20) change in probability.)
								int g2 = g;
								if (!g) g2 = (1 << 15) - 1;

								int oldshift = shiftflagmode;
								rawvals[g][oldshift] = exp(dous[j]->doanalyze<aroundturner>(tb, turn, chromstarts[i],
									chromstarts[i + 1] - 1, q, -1, -1, true, 0, -20 + factor) - factor);
								shiftflagmode = oldshift;

								if (rawvals[g][oldshift] < 0) continue;

								for (int t = 0; t < TYPEBITS + 1; t++)
								{
									if (g2 & (1 << t))
									{
										sumnegval[t] += rawvals[g][shiftflagmode];
									}
								}

							}
						}

						for (int t = 0; t < TYPEBITS + 1; t++)
						{
							if (!sumnegval[t]) sumnegval[t] = 1;
						}


#pragma omp critical(negshifts)
						{
							for (int g = 0; g < NUMTYPES * 2; g++)
							{
								for (int s = 0; s < NUMSHIFTS; s++)
								{
									double val = rawvals[g][s];

									if (_finite(val) && val > 1e-10)
									{
										int g2 = g;
										if (!g) g2 = (1 << 15) - 1;							      

										// This is hardcoded for the generation count of 3.
										int marker = -q - 1000;
										dous[j]->negshift[marker] += val * (1.0 - ((g >> 6) & 1) * 2) * ((g2 >> 6) & 1) / sumnegval[6];

										if (dous[j]->pars[0])
											dous[j]->pars[0]->negshift[marker] += val * (1.0 - ((g >> 0) & 1) * 2) * ((g2 >> 0) & 1) / sumnegval[0];

										if (dous[j]->pars[1])
											dous[j]->pars[1]->negshift[marker] += val * (1.0 - ((g >> 3) & 1) * 2) * ((g2 >> 3) & 1) / sumnegval[3];
										if (dous[j]->gen >= 2)
										{
											if (dous[j]->pars[0] && dous[j]->pars[0]->pars[0])
												dous[j]->pars[0]->pars[0]->negshift[marker] += val * (1.0 - ((g >> 1) & 1) * 2) * ((g2 >> 1) & 1) / sumnegval[1];

											if (dous[j]->pars[0] && dous[j]->pars[0]->pars[1])
												dous[j]->pars[0]->pars[1]->negshift[marker] += val * (1.0 - ((g >> 2) & 1) * 2) * ((g2 >> 2) & 1) / sumnegval[2];

											if (dous[j]->pars[1] && dous[j]->pars[1]->pars[0])
												dous[j]->pars[1]->pars[0]->negshift[marker] += val * (1.0 - ((g >> 4) & 1) * 2) * ((g2 >> 4) & 1) / sumnegval[4];

											if (dous[j]->pars[1] && dous[j]->pars[1]->pars[1])
												dous[j]->pars[1]->pars[1]->negshift[marker] += val * (1.0 - ((g >> 5) & 1) * 2) * ((g2 >> 5) & 1) / sumnegval[5];
										}
									}
								}
							}
						}
					}					

					// Print the results for the four cases, if we keep track of them
					if (doprint)
					{
						double sum = 0;
						if (skipsome)
						{
							
							if (fabs(probs[0] + probs[1] - probs[2] - probs[3]) <
								fabs(probs[0] + probs[2] - probs[1] - probs[3]))
							{
								probs[1] = 0;
							}
							else 
							{
								probs[2] = 0;
							}
							if (probs[0] > probs[3])
							{
								probs[3] = 0;
							}
							else
							{
								probs[0] = 0;
							}
						    /*						    const double IMPOSSVAL = 2;
						    for (int j = 0; j < 2; j++)
						      {
							double min = IMPOSSVAL;
							int mini = 0;
							for (int i = 0; i < 4; i++)
							  {
							    if (probs[i] < min)
							      {
								min = probs[i];
								mini = i;
							      }
							  }

							probs[mini] = IMPOSSVAL;
						      }

						    for (int i = 0; i < 4; i++)
						      {
							if (probs[i] == IMPOSSVAL) probs[i] = 0;
							}*/
						  }

						for (int i = 0; i < 4; i++)
						{
							sum += probs[i];
						}

						if (sum == 0)
						{
							sum = 1;
						}
						for (int i = 0; i < 4; i++)
						{
							probs[i] /= sum;
						}

						sprintf(lineout, "%lf\t%lf\t%lf\t%lf\n", probs[0], probs[1], probs[2], probs[3]);
						//sprintf(lineout, "HEHE R");
						if (oqp[j] > 50000) printf("%d\t%d\n", j, oqp[j]);
						strcpy(&outqueue[j][oqp[j]], lineout);
						oqp[j] += strlen(lineout);
					}

					if (!full)
					{
						int marker = -q - 1000;

						// Contribute haplotype data, but truncate it, i.e. a 50/50 contribution for either interpretation is not added.
						// Instead, we have a cap later on at the maximum change at any iteration.
						{
#pragma ivdep
							for (int k = 0; k < (int) reltree.size(); k++)
							{
								int i = reltree[k]->n;
								if (haplos[i][0] || haplos[i][1])
								{
									float base;
									if (reltree[k] != dous[j])
									{
										base = min(haplos[i][0], haplos[i][1]);
									}
									else
										base = 0;
#pragma omp critical(update)
									{
										getind(i)->haplobase[marker] += haplos[i][0] - base;
										getind(i)->haplocount[marker] += haplos[i][1] + haplos[i][0] - base * 2;
									}
									haplos[i][0] = 0;
									haplos[i][1] = 0;
								}
							}
						}
					}
				}
			}
		}

		printf("Printing output\n");
		for (unsigned int j = 0; j < outqueue.size(); j++)
		{
			fprintf(out, "%d\n", dous[j]->n);
			fprintf(out, "%s\n", &outqueue[j].front());
		}
	}

	if (!full)
	{
		for (unsigned int i = 0; i < 2000000; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;

			for (unsigned int j = 0; j < ind->haplocount.size(); j++)
			{
				if (ind->haplocount[j] && ind->haploweight[j])
				{
					double b1 = ind->haplobase[j];
					double b2 = ind->haplocount[j] - ind->haplobase[j];

					b1 /= ind->haploweight[j];
					b2 /= (1.0 - ind->haploweight[j]);
					double intended = b1 / (b1 + b2);

					//printf("---\n");
					/*for (int q = 0; q < 20; q++)
					{
					double b12 = b1 * intended;
					double b22 = b2 * (1.0 - intended);
					intended = b12 / (b12 + b22);
					printf("%lf\n", intended);
					}*/

					//				double intended = ind->haplobase[j] / ind->haplocount[j];
					//double intended2 = intended * (1.0 - ind->haploweight[j]) +
					//(1.0 - intended) * ind->haploweight[j];

					//if (early)
					{
						double nnn = 3;
						double limn = (nnn - 1.0) * ind->haploweight[j] * (-1 + ind->haploweight[j]);

						double limd1 = -1 - (nnn - 1.0) * ind->haploweight[j];
						double limd2 = (nnn - 1.0) * ind->haploweight[j] - nnn;

						double lim = min(limn/limd1, limn/limd2);

						double diff = intended - ind->haploweight[j];

						if (fabs(diff) > lim)
						{
							intended = ind->haploweight[j] + diff / fabs(diff) * lim;
						}
					}
					intended = min(intended, 0.999);
					ind->haploweight[j] = max(intended, 0.001);
				}
			}

			if (ind->pars[0] || ind->pars[1] || !ind->haplocount.size()) continue;

			// Perform the inversions indicated by the negshift data, at most a single one per individual
			// and chromosome, maybe 0.
			for (int c = 0; c < (int) chromstarts.size() - 1; c++)
			{
				int minstart = chromstarts[c + 1];
				double minval = -0.2;

				for (int p = chromstarts[c]; p < (int) chromstarts[c + 1]; p++)
				{
					if (ind->negshift[p] < minval)
					  {
						  minstart = p;
						  minval = ind->negshift[p];
					  }
				}

				for (int p = minstart + 1; p < (int) chromstarts[c + 1]; p++)
				{
					ind->haploweight[p] = 1.0f - ind->haploweight[p];
				}
			}
		}
	}
	generation++;
}

// Handle the fuss of trailing \r, \n characters when combining scanf and gets and stuff.
void clean(char* tlf)
{
	int i = strlen(tlf);
	while (i && tlf[i - 1] < 32)
	{
		tlf[i - 1] = 0;
		i--;
	}
}

int main()
{
#ifdef _MSC_VER
	// Performant?
	// This turns off some aspects of IEEE precision, but we do not really need it, as we avoid denormalized values
	// anyway.
	_controlfp(_DN_FLUSH, _MCW_DN);
#endif

	printf("Discretization step: ");
	scanf("%lf", &discstep);

	printf("Number of sexes in map: ");	
	scanf("%d", &sexc);

	// Not really related to the number of generations, but doing it like this makes the
	// estimates similar for the two cases. Whether the old 4-state, 2-gen model was
	// actually correct is another issue entirely.
	//
	// / 50 everywhere is probably more appropriate
	if (NUMGEN == 3)
	{
		baserec[0] = -discstep / 50.0;
	}
	else
	{
		baserec[0] = -discstep / 50.0;
	}
	baserec[1] = baserec[0];

	char tlf[255];
	fgets(tlf, 255, stdin);

	printf("Marker info file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);
	FILE* in = fopen(tlf, "r");
	readmarkerinfo(in);
	fclose(in);

	printf("Pedigree file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);
	in = fopen(tlf, "r");
	readped(in);
	fclose(in);

	printf("Marker data file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);
	in = fopen(tlf, "r");
	readmarkerdata(in);
	fclose(in);

	/*readqtlmas();*/

	//sort(dous.begin(), dous.end());

	printf("Output file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);

	postmarkerdata();

	FILE* out = fopen(tlf, "w");

	if (HAPLOTYPING)
	for (int i = 0; i < 23; i++)
	{
		//		  	  	{
		early = (i < 1);
		doit<false>(out);
		//		}

		//	doit<true>(out);
		for (unsigned int i = 0; i < 2000000; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;

			if (ind->haplocount.size())
			{
				fprintf(out, "%d\n", i);
				// Printing of haplotype data for each iteration
				for (unsigned int j = 0; j < chromstarts[1]; j++)
				{
					fprintf(out, "%f\t%d\t%d\t\t%f\n", ind->haploweight[j], ind->markerdata[j].first, ind->markerdata[j].second, ind->negshift[j]);
					ind->negshift[j] = 0;
				}
				fprintf(out, "\n");
			}
		}
	}
	doit<true>(out);
	fclose(out);	
}
