// Master sorting v6 for Approximate Sorting experiment
// Author: Aditya Narasimhan


// #include <bits/stdc++.h>
#include <string>
#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "err.h"
#include <signal.h>
#include <chrono> //for sleep time
#include <thread>
#include <cstdlib> //for abs
#include <vector> //for generating a permutation of numElements
#include <math.h>       // for the ceil function

using namespace std;

// GLobal elements
unsigned long int numElements = 0;
unsigned long int* arr; //array of length numElements
unsigned long int* y;
unsigned long int* yIndex;
unsigned long int* tempArr = new unsigned long int[5];
unsigned long int* tempArrIndex = new unsigned long int[5];
unsigned long int pivotIndexActual = 0;
long int maxDifference = 0;
long int recursiveCount=0;
long int maxD = 0;


inline void swap(unsigned long int* a, unsigned long int* b) {
	unsigned long int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

// Function to generate n non-repeating random numbers 
void generateRandom() 
{
	for(long int i = 0; i < numElements; i++)
		// arr[i] = int(rand()) % 1000;
		arr[i] = int(rand()) % 10000000000000;
} 

/* A utility function to print array of size n */
void printArray() 
{ 
	for (long int i=0; i<numElements; ++i) 
		cout << arr[i] << " "; 
	cout << "\n"; 
} 
 
/* A utility function to print array of size n without giving n */
void printArray_v2(unsigned long int* arrTemp) 
{ 
	for (long int i=0; i<numElements; ++i) 
		cout << i << ":" << arrTemp[i] << " "; 
	cout << "\n"; 
} 

/* A utility function to print array of size n without giving n */
void printArray_v3(unsigned long int* arrTemp, long int l) 
{ 
	for (long int i=0; i < l; ++i) 
		cout << arrTemp[i] << " "; 
	cout << "\n"; 
} 


inline unsigned long int median_1(unsigned long int a[], unsigned long int aInd[]) {
	pivotIndexActual = aInd[0];
	return a[0];
}

inline unsigned long int median_2(unsigned long int a[], unsigned long int aInd[]) {
	pivotIndexActual = aInd[1];
	return a[1];
}

inline unsigned long int median_3(unsigned long int a[], unsigned long int aInd[]) {
	unsigned long int temp;
	if (a[0] < a[1]) {
		swap(&a[0], &a[1]);
		swap(&aInd[0], &aInd[1]);
	}
	if (a[0] > a[2]) {
		a[0] = a[1];
		a[1] = a[2];

		aInd[0] = aInd[1];
		aInd[1] = aInd[2];
	} else {
		pivotIndexActual = aInd[0];
		return a[0];
	}
	if (a[0] > a[1]) {
		pivotIndexActual = aInd[0];
		return a[0];
	} else {
		pivotIndexActual = aInd[1];
		return a[1];
	}
}

inline unsigned long int median_4(unsigned long int a[], unsigned long int aInd[]) {
	unsigned long int temp1, temp2, temp1Ind, temp2Ind;
	if (a[0] < a[1]) {
		swap(&a[0], &a[1]);
		swap(&aInd[0], &aInd[1]);
	}
	if (a[2] > a[3]) {
		swap(&a[1], &a[2]);
		swap(&aInd[1], &aInd[2]);
	} else {
		temp1 = a[1];
		temp2 = a[2];
		a[1] = a[3];
		a[2] = temp1;
		a[3] = temp2;

		temp1Ind = aInd[1];
		temp2Ind = aInd[2];
		aInd[1] = aInd[3];
		aInd[2] = temp1Ind;
		aInd[3] = temp2Ind;
	}
	if (a[0] > a[1]) {
		a[0] = a[1];
		a[1] = a[2];

		aInd[0] = aInd[1];
		aInd[1] = aInd[2];
	} else {
		a[1] = a[3];

		aInd[1] = aInd[3];
	}
	if (a[0] > a[1]) {
		pivotIndexActual = aInd[0];
		return a[0];
	} else {
		pivotIndexActual = aInd[1];
		return a[1];
	}
}

inline unsigned long int median_5(unsigned long int a[], unsigned long int aInd[]) {
	unsigned long int temp1, temp2, temp1Ind, temp2Ind;
	if (a[0] < a[1]) {
		swap(&a[0], &a[1]);
		swap(&aInd[0], &aInd[1]);
	} 
	if (a[2] > a[3]) {
		swap(&a[1], &a[2]);
		swap(&aInd[1], &aInd[2]);
	} else {
		temp1 = a[1];
		temp2 = a[2];
		a[1] = a[3];
		a[2] = temp1;
		a[3] = temp2;

		temp1Ind = aInd[1];
		temp2Ind = aInd[2];
		aInd[1] = aInd[3];
		aInd[2] = temp1Ind;
		aInd[3] = temp2Ind;
	}
	if (a[0] > a[1]) {
		a[0] = a[1];
		a[1] = a[3];
		a[3] = a[4];

		aInd[0] = aInd[1];
		aInd[1] = aInd[3];
		aInd[3] = aInd[4];
	} else {
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];

		aInd[1] = aInd[2];
		aInd[2] = aInd[3];
		aInd[3] = aInd[4];		
	}
	if (a[2] > a[3]) {
		swap(&a[1], &a[2]);
		swap(&aInd[1], &aInd[2]);
	} else {
		temp1 = a[1];
		temp2 = a[2];
		a[1] = a[3];
		a[2] = temp1;
		a[3] = temp2;

		temp1Ind = aInd[1];
		temp2Ind = aInd[2];
		aInd[1] = aInd[3];
		aInd[2] = temp1Ind;
		aInd[3] = temp2Ind;
	}
	if (a[0] > a[1]) {
		a[0] = a[1];
		a[1] = a[2];

		aInd[0] = aInd[1];
		aInd[1] = aInd[2];
	} else {
		a[1] = a[3]; 

		aInd[1] = aInd[3]; 
	}
	if(a[0] > a[1]) {
		pivotIndexActual = aInd[0];
		return a[0];
	} else {
		pivotIndexActual = aInd[1];
		return a[1];
	}
}

inline unsigned long int median(unsigned long int aTemp[], unsigned long int aTempIndex[], long int n)
{
	switch(n) 
    {
		case 1: 
			return median_1(aTemp, aTempIndex);
			break;
		case 2:
			return median_2(aTemp, aTempIndex);
			break;
		case 3:
			return median_3(aTemp, aTempIndex);
			break;
		case 4:
			return median_4(aTemp, aTempIndex);
			break;
		case 5:
			return median_5(aTemp, aTempIndex);
			break;
		default:
			break;
	}
	return -1;
}


void copyLessThan5(unsigned long int* aPointer, long int startInd, long int num)
{
	// cout << "Subarray: ";
	for(long int i = 0; i < num; i++)
	{
		// cout << "[" << startInd + i << "]: ";
		// cout << aPointer[startInd + i] << " ";
		tempArr[i] = aPointer[startInd + i];
		tempArrIndex[i] = startInd + i;
	} // for i
	// cout << endl;
}

void copyLessThan5yIndex(unsigned long int* aPointer, long int startInd, long int num)
{
	// cout << "Subarray Y: ";
	for(long int i = 0; i < num; i++)
	{
		// cout << "[" << yIndex[startInd + i] << "]: ";
		// cout << aPointer[startInd + i] << " ";
		tempArr[i] = aPointer[startInd + i];
		tempArrIndex[i] = yIndex[startInd + i];
	} // for i
	// cout << endl;
}

void findMedian(long int n1, long int n2, long int lIndex, long int rIndex)
{
	long int t;
	long int pos;
	unsigned long int tempSum;

	for (int i=0; i < n1; i++) 
	{
		y[i] = 0;
        // y[i] = y[i] + findSum(x, lIndex+i*5, 5);
		// cout << "i: " << i << endl;
		copyLessThan5(arr, lIndex+i*5, 5);
		y[i] = y[i] + median(tempArr, tempArrIndex, 5);
		yIndex[i] = pivotIndexActual;
		// cout << " Median: " << y[i] << " in position " << yIndex[i] << endl;
		// cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}

	y[n1] = 0;
	copyLessThan5(arr, lIndex+n1*5, n2);
    y[n1] = y[n1] + median(tempArr, tempArrIndex, n2);
	// y[n1] = y[n1] + findSum(x, lIndex+(n1*5), n2);
	yIndex[n1] = pivotIndexActual;
	// cout << " Median: " << y[n1] << " in position " << yIndex[n1] << endl;

    if (n2 > 0) t = n1+1;
    else t = n1;

    // cout << "n1 = " << n1 << "; n2 = " << n2 << endl;
    // for (int i=0; i < t; i++)
    // 	cout << y[i] << ", ";
    // cout << endl;
    // cout << "-----------------" << endl;

	while (t > 1) {
		n1 = t/5;
		n2 = t % 5;
		if (n1 == 0) {
			pos = 0;
			tempSum = 0;
			// for (int i=0; i < n2; i++)
			// 	tempSum = tempSum + y[i]; //y[i] should have the median
			copyLessThan5yIndex(y, 0, n2);
			tempSum = tempSum + median(tempArr, tempArrIndex, n2);
			y[pos] = tempSum;
			yIndex[pos] = pivotIndexActual;
		}
		else {
			for (int i=0; i < n1; i++) {
				tempSum = 0;
                // tempSum = tempSum + findSum(y, i*5, 5);
				copyLessThan5yIndex(y, i*5, 5);
				tempSum = tempSum + median(tempArr, tempArrIndex, 5);
				y[i] = tempSum;
				yIndex[i] = pivotIndexActual;
			}
			tempSum = 0;
            // tempSum = tempSum + findSum(y, n1*5, 5);
			copyLessThan5yIndex(y, n1*5, n2);
			tempSum = tempSum + median(tempArr, tempArrIndex, n2);
			y[n1] = tempSum;
			yIndex[n1] = pivotIndexActual;
		}

		if (n2 > 0) t = n1+1;
    	else t = n1;

    	// cout << "n1 = " << n1 << "; n2 = " << n2 << endl;
	    // for (int i=0; i < t; i++)
	    // 	cout << y[i] << ", ";
	    // cout << endl;
	    // cout << "-----------------" << endl;
	} // while

}


// used to find the number of inversions
long int numInversions()
{
    long int countInversions = 0;
    for(long int i = 0; i < numElements; i++)
    {
        for(long int j = i + 1; j < numElements; j++) // looping through all elements to the right
        {
            if(arr[j] < arr[i]) // checking for the number of elements lesser than the current element to the right
                countInversions++;
        }
    }
    return countInversions;
}

// FULL BINARY SEARCH used for searching for element in given arr for MD calculating
long int binarySearch(unsigned long int* arr, long int l, long int r, long int x)
{
    if (r >= l) {
        int mid = l + (r - l) / 2;

        // If the element is present at the middle
        // itself
        if (arr[mid] == x)
            return mid;

        // If element is smaller than mid, then
        // it can only be present in left subarray
        if (arr[mid] > x)
            return binarySearch(arr, l, mid - 1, x);

        // Else the element can only be present
        // in right subarray
        return binarySearch(arr, mid + 1, r, x);
    }

    // We reach here when element is not
    // present in array
    return -1;
}

// used to find the maximum displacement or the Chebyshev distance
long int chebyshevDistance(unsigned long int* sortedArr)
{
    long int maxDiff = 0;
    for(long int i = 0; i < numElements; i++)
    {
        long int originalPos = binarySearch(sortedArr, 0, numElements, arr[i]); // search for the element arr[i] in
                                                                // the sorted array and returns the rank/original position
        if(originalPos == -1)
        {
            cout << "Element not found" << endl;
            return -1;
        }
        // cout << "away: " << arr[i] << " is in " << i << " but has to be in " << originalPos <<" => Diff = " << abs(originalPos - i) << endl;

        // check for max
        int currDist = originalPos - i;
        if(currDist < 0) {
			currDist = currDist + (currDist * 2); //take the absolute value
		}
        if(currDist > maxDiff)
            maxDiff = abs(originalPos - i);
    }

    return maxDiff;
}

// FULL partition algorithm for quick sort
long int partition (unsigned long int* arr, long int low, long int high) 
{ 
    long int pivot = arr[high]; // pivot 
    long int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
  
    for (long int j = low; j <= high - 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[j] < pivot) 
        {
            i++; // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        }
    }
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 

// FULL quick sort algorithm
void quickSort(unsigned long int* arr, long int low, long int high) 
{
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
        at right place */
        long int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 

// Partition for MedianRearrange
long int Partition(unsigned long int* a, long int low, long int high)
{
	long int pivot, index, i;
	index = low;
	pivot = high;
  
	// Getting index of the pivot.
	for(long int i=low; i < high; i++)
	{
		if(a[i] < a[pivot])
		{
			swap(&a[i], &a[index]);
			index++;
		}
	}
	// Swapping value at high and at the index obtained.
	swap(&a[pivot], &a[index]);
 
	return index;
}

// MEDIAN selection of pivot for MedianRearrange
long int RandomPivotPartition(unsigned long int* a, long int low, long int high)
{
	long int n, temp;
	long int n1, n2;
	// n = rand();
	// Randomizing the pivot value in the given subpart of array.
	int l = low, r = high; // index positions as limits inclusive
	n1 = (r-l+1)/5;
	n2 = (r-l+1)%5;
	
	if (n2 > 0) 
	{
		yIndex = new unsigned long int [n1+1];
		y = new unsigned long int [n1+1];
	}
	else 
	{
		yIndex = new unsigned long int [n1+1];
		y = new unsigned long int [n1+1];
	}

    findMedian(n1, n2, l, r); // next step is to add index - DONE
    unsigned long int pivot = y[0];
	unsigned long int pivotIndex = yIndex[0];
	// cout << "MEDIAN: " << y[0] << endl;
    // to find the index of pivot

	// cout << "COMPARE - Median Pivot is arr[" << pivotIndex << "]: " << pivot << endl;
 
	// Swapping pivot value from high, so pivot value will be taken as a pivot while partitioning.
	swap(&a[high], &a[pivotIndex]);
 
	return Partition(a, low, high);
}

// MedianRearrange

long int MedianRearrangeApproxSort(unsigned long int* a, long int low, long int high, long int D)
{
	long int pindex;
	if(low < high && D <= maxD)
	{
		// Partitioning array using randomized pivot.
		pindex = RandomPivotPartition(a, low, high);

		// cout << "After partition: ";
		// printArray();
        // cout << "D: " << D << endl;
		
		// Recursively implementing MedianRearrangeApproxSort.
		MedianRearrangeApproxSort(a, low, pindex-1, D+1);
        // cout << "D inbetween: " << D << endl;
		MedianRearrangeApproxSort(a, pindex+1, high, ++D);
        // cout << "D out: " << D << endl;
		// cout << "Number of inversions: " << numInversions() << endl;
		// D++;
	}

	return 0;
}

bool isSorted(unsigned long int* array, long int len)
{
	if(len < 2)
	{
		// cout << "returning 0" << endl;
		return true;
	}
	
	int previous = array[0];
	
	for(long int i = 1; i < len; i++)
	{

		if(previous <= array[i])
		{
			previous = array[i];
		}
		else
		{
			// cout << "returning 0" << endl;
			// cout << "NOT SORTED ARRAY: ";
			// printArray_v3(array, len);
			return false;

		}
	}
	
	// cout << "returning 1" << endl;
	// cout << "SORTED ARRAY: ";
	// printArray_v3(array, len);
	return true;
}

unsigned long int* mergeArrays(unsigned long int* arr1, unsigned long int* arr2, long int n1, long int n2)
{
	unsigned long int* arr3 = new unsigned long int[n1+n2];

	// for(long int i = 0; i < n1; i++)
	// 	cout << "arr1["<<i<<"] = "<<arr1[i]<<endl;
	
	// cout << endl;

	// for(long int i = 0; i < n2; i++)
	// 	cout << "arr2["<<i<<"] = "<<arr2[i]<<endl;

	// cout << endl;

	long int i = 0, j = 0, k = 0; 
  
    // Traverse both array 
    while (i<n1 && j <n2) 
    { 
        // Check if current element of first 
        // array is smaller than current element 
        // of second array. If yes, store first 
        // array element and increment first array 
        // index. Otherwise do same with second array 
        if (arr1[i] < arr2[j]) 
            arr3[k++] = arr1[i++]; 
        else
            arr3[k++] = arr2[j++]; 
    }
    // Store remaining elements of first array 
    while (i < n1) 
        arr3[k++] = arr1[i++]; 
  
    // Store remaining elements of second array 
    while (j < n2) 
        arr3[k++] = arr2[j++]; 

	return arr3;
}

unsigned long int* adaptiveSort(unsigned long int* a, long int length)
{
	recursiveCount++;
	if(length < 2) return a;

	//split to odd and even arrays 
	// cout << "CEIl LENGTH: " << int(ceil(length/2.0)) << endl;
	// cout << "floor LENGTH: " << int(floor(length/2.0)) << endl;
	unsigned long int* evens = new unsigned long int[int(ceil(length/2.0))];
	unsigned long int* odds = new unsigned long int[int(floor(length/2.0))];

	int lim = length / 2;

	for(long int i = 0, i2 = 0; i < lim; i++, i2++)
	{
		evens[i] = a[i2];
		i2++;
		odds[i] = a[i2];
	}

	if(length % 2 > 0)
	{
		evens[(length/2)] = a[length - 1];
		//cout << "evens[(length/2)]: " << evens[(length/2)] << endl;
	}

	// cout << "evens\n";
	// printArray_v3(evens, int(ceil(length/2.0)));
	// cout << "odds\n";
	// printArray_v3(odds, int(floor(length/2.0)));

	if(!isSorted(evens, int(ceil(length/2.0))))
	{
		evens = adaptiveSort(evens, int(ceil(length/2.0)));
	}
	if(!isSorted(odds, int(floor(length/2.0))))
	{
		// cout << "NOT SORTED odds: ";
		// printArray_v3(odds, int(floor(length/2.0)));
		odds = adaptiveSort(odds, int(floor(length/2.0)));
	}

	// cout << "evens\n";
	// printArray_v3(evens, int(ceil(length/2)));
	// cout << "odds\n";
	// printArray_v3(odds, int(ceil(length/2)));

	a = mergeArrays(evens, odds, int(ceil(length/2.0)), int(floor(length/2.0)));

	// cout << "Completion: Sorted array:\n";	
	// printArray_v3(a, int(ceil(length)));

	delete [] evens;
	delete [] odds;

	return a;
}

// used for searching for element in given arr for MD
long int modbinarySearch(unsigned long int* arr, long int l, long int r, long int x, long int md)
{
	if((r - l + 1) <= (4 * md)) // length becomes less than 4L
	{
		for(long int z = l; z <= r; z++)
		{
			if(arr[z] == x)
				return z; // return index position if found
		}
		return -1; // not found, so return -1
	}

    if (r >= l) 
	{
        long int mid = l + (r - l) / 2;

        // If the element is present at the middle
        // itself
        if (arr[mid] == x)
            return mid;

        // If element is smaller than mid, then
        // it can only be present in left subarray
        if (arr[mid] > x)
            return modbinarySearch(arr, l, mid + (2*md)  - 1, x, md);

        // Else the element can only be present
        // in right subarray
        return modbinarySearch(arr, mid - (2*md) + 1, r, x, md);
    }

    // We reach here when element is not
    // present in array
    return -1;
}

//main function
int main(int argc, char **argv)
{   
	numElements = atoi(argv[1]);
	long int numTimes = atoi(argv[2]);
	maxD = atoi(argv[3]);

	arr = new unsigned long int[numElements];

	// unsigned long int* sortedArray;
	// unsigned long int* sortedArrayMD;
	clock_t start1, end1;
	clock_t start2, end2;
	clock_t start3, end3;
	double* timeAvg1 = new double[numTimes];
	double* timeAvg2 = new double[numTimes];
	double* timeAvgSearch = new double[10000];

	cout << "NUMELEMENTS: " << numElements << endl;

	for(long int avgTimes = 0; avgTimes < numTimes; avgTimes++)
	{
		// cout << "iteration: " << avgTimes << endl;
		generateRandom(); 

		// cout << "Median in index position " << numElements/2 << endl;
		// cout << "Median is " << select((numElements/2) + 1, arr, numElements) << endl;

		// cout << "Randomly generated array " << avgTimes << endl;
		// printArray();

		// adaptive sort
		start1 = clock();
		// sortedArray = adaptiveSort(arr, numElements);
		end1 = clock();
		double time_taken1 = double(end1 - start1)/CLOCKS_PER_SEC;
		timeAvg1[avgTimes] = time_taken1;

		// cout << fixed << timeAvg1[avgTimes] << setprecision(5) << " - A" << endl;
		
		// cout << "Recursive count: " << recursiveCount << endl;

		// cout << "SORTED randomly generated array: " << endl;
		// printArray_v2(sortedArray);

		// unsigned long int randomCheckNum = 1 + ( std::rand() % ( numElements - 1 + 1 ) );
		// cout << "Adaptive sort: Random check if sorted: arr[" << randomCheckNum << "]: " << sortedArray[randomCheckNum] << endl;

		// quick sort
		start2 = clock();
		// quickSort(arr, 0, numElements-1); // full quick sort
		MedianRearrangeApproxSort(arr, 0, numElements-1, 0); // median rearrange
		end2 = clock();
		double time_taken2 = double(end2 - start2)/CLOCKS_PER_SEC;
		timeAvg2[avgTimes] = time_taken2;

		// cout << fixed << timeAvg2[avgTimes] << setprecision(5) << " - Q" << endl;
		
		// cout << "Quick sort: Random check if sorted: arr[" << randomCheckNum << "]: " << arr[randomCheckNum] << endl;

	}

	double total1 = 0.0, total2 = 0.0, total3 = 0.0;
	for(long int avgTimes = 0; avgTimes < numTimes; avgTimes++)
	{
		total1 += timeAvg1[avgTimes];
		total2 += timeAvg2[avgTimes];		
	}

	// cout << "Average adaptive sort time: " << total1/numTimes << endl;
	cout << "Average quick sort time: " << total2/numTimes << endl;

	// cout << "Number of inversions: " << numInversions() << endl;
	
	// cout << "End of approximately sorting: " << endl;
	// printArray_v2(arr);

	// sortedArrayMD = adaptiveSort(arr, numElements);
	// cout << "FULLY SORTED ARRAY: " << endl;
	// printArray_v2(sortedArrayMD);
	
	// int maxDisp = chebyshevDistance(sortedArrayMD);
	long int maxDisp = ceil((numElements)/pow(2,(maxD+1)));
	cout << "Maximum displacement ceiling: " << maxDisp << endl;

	for(long int b = 0; b < 10000; b++)
	{
		unsigned long int toSearch = arr[int(rand()) % (numElements - 1)];
		start3 = clock();
		long int ifFound = modbinarySearch(arr, 0, numElements - 1, toSearch, maxDisp);
        // long int ifFound = binarySearch(arr, 0, numElements, toSearch);
		end3 = clock();
		double time_taken3 = double(end3 - start3)/CLOCKS_PER_SEC;
		timeAvgSearch[b] = time_taken3;
		// if(ifFound != -1)
		// 	cout << "Element " << toSearch << " found at " << ifFound << endl;
		// else
		// 	cout << "Element " << toSearch << " NOT found" << endl;
	}

	for(long int b = 0; b < 10000; b++)
	{
		total3 += timeAvgSearch[b];	
	}

	cout << "Average modified search time for " << 10000 << " elements: " << total3/10000 << endl;

    delete [] arr;
	// delete [] sortedArrayMD;
	// delete [] sortedArray;
    delete [] y;
    delete [] yIndex;

	return 0;
}

// Approximately sorted array: 
// 1 3 5 4 6 8 9 10 11 7 12 13 14 15 2 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
// The approximate array's maximum fake rank is: 13

// what is the point of implementing a O(n) algorithm? 
// Median rearrange is still recursive
// How about submitting to SWAT itself? 
