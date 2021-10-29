classdef HeapBasedPriorityQueue
    
    
    
    %   arrayBuild = [23, 34, 45, 45, 1, 33, 99, 78, 786, 546,333, 423, 876]';
    %   heapSize = length(arrayBuild);
    
    %  arrayBuild = buildMaxHeap(arrayBuild, heapSize);
    
    %   [arrayBuild, heapSize]  = insert(arrayBuild, 23, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 3, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 76, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 1298, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 12, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 22, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 2244, heapSize);
    %   [arrayBuild, heapSize] = insert(arrayBuild, 0222, heapSize);
    
    
    methods(Static)
      
        function [builtArray, keepIndexes] = buildMaxBasedHeapOnArray(arrayToBuild, heapSize)
            arrayIndexes(1:length(arrayToBuild)) = 1:length(arrayToBuild);
            [builtArray,  keepIndexes, ~] = HeapBasedPriorityQueue.buildMaxHeap(arrayToBuild, arrayIndexes, heapSize);
        end
        
        function [arrayBuild, arrayIndexes] = insertEleInArray(arrayToBuild, ele, heapSize)
            arrayIndexes(1:length(arrayToBuild)) = 1:length(arrayToBuild);
            [arrayBuild, ~, arrayIndexes]  = HeapBasedPriorityQueue.insert(arrayToBuild, ele, heapSize, arrayIndexes);
            arrayBuild = arrayBuild(1:heapSize);
            arrayIndexes = arrayIndexes(1:heapSize);
        end
        
        function retVal = getRightChild(array, index)
            retVal = -1;
            if ( (((2*index)+1) <= length(array))  &&  (index >= 2) )  % >= 2 because in matlab the index starts from 1 unlike 0 in C
                retVal = ((2*index)+1);
                return;
            end
            return;
        end
        
        function retVal = getLeftChild(array, index)
            retVal = -1;
            if ( ((2*index) <= length(array))  &&  (index >= 2) )  % >= 2 because in matlab the index starts from 1 unlike 0 in C
                retVal = (2*index);
                return;
            end
            return;
        end
        
        function retVal = getParent(array, index)
            retVal = -1;
            if ( ( index <= length(array))  &&  (index > 2) )  % >= 2 because in matlab the index starts from 1 unlike 0 in C
                retVal = floor(index/2);
                return;
            end
            return;
        end
        
        function [array, keepIndexes, heapSize] = maxHeapify(array, index, keepIndexes, heapSize)
            leftChildIndex = 2*index;   % getLeftChild(array, index);
            rightChildIndex = ((2*index)+1); % getRightChild(array, index);
            
            % finding largest among index, left child and right child
            largest = index;
            
            if ( (leftChildIndex <= heapSize) && (leftChildIndex > 1) && ...
                    (array(leftChildIndex) > array(largest))) % for max heap the root will always contain bigger value than the child nodes
                largest = leftChildIndex;
            end
            
            if ( ((rightChildIndex <= heapSize) && (rightChildIndex > 1)) && ....
                    (array(rightChildIndex) > array(largest)) )  % for max heap the root will always contain bigger value than the child nodes
                largest = rightChildIndex;
            end
            
            
            % largest is not the node, node is not a heap
            if (largest ~= index)
                % swapping
                temp = array(largest);
                array(largest) = array(index);
                array(index) = temp;

                tempIndex = keepIndexes(largest);
                keepIndexes(largest) = keepIndexes(index);
                keepIndexes(index) = tempIndex;

                [array, keepIndexes, heapSize] = HeapBasedPriorityQueue.maxHeapify(array, largest, keepIndexes, heapSize);
            end
            return;
        end
        
        function [array,  keepIndexes, heapSize] = buildMaxHeap(array, keepIndexes, heapSize)
            for i= floor(heapSize/2):-1:1
                [array, keepIndexes, heapSize] = HeapBasedPriorityQueue.maxHeapify(array, i, keepIndexes, heapSize);
            end
            return;
        end
        
        function retVal = getMaximumVals(array)
            retVal = array(2);
            return;
        end
        
        function [maxm, heapSize] = extractMax(array, heapSize)
            maxm = array(2);
            array(2) = array(heapSize); % putting the last element of the array at the first location
            heapSize = hipSize-1;
            max_heapify(array, 1);
            return;
        end
        
        function [array, keepIndexes] = increaseKey(array, index, key, keepIndexes)
            array(index) = key;
            keepIndexes(index) = index;
            while( (index > 1) && (array(floor(index/2)) < array(index)) ) % the reason to check here (index > 2) is that if there are only 2 elmenet then the biggest element will be at the first cell and the other element will be at the second cell. So, there is nothing to making heap sort
                temp = array(floor(index/2));
                array(floor(index/2)) = array(index);
                array(index) = temp;
                
                tempIndex = keepIndexes(floor(index/2));
                keepIndexes(floor(index/2)) = keepIndexes(index);
                keepIndexes(index) = tempIndex;
                
                index = floor(index/2);
            end
            return;
        end
        
        function decreaseKey(array, index, key)
            array(index) = key;
            HeapBasedPriorityQueue.maxHeapify(array, index);
        end
        
        function [array, heapSize, keepIndexes] = insert(array, key, heapSize, keepIndexes)
            INF = 100000;
            heapSize = heapSize + 1;
            array(heapSize) = (-1 * INF);
            [array, keepIndexes] = HeapBasedPriorityQueue.increaseKey(array, heapSize, key, keepIndexes);
            return;
        end
    end
end
