## Sat 09/16/2023
- bit arithmetic approach to swap
    ```cpp
    void swap(int* a, int* b) {
        a ^= b;
        b ^= a;
        a ^= b;
    }
    ```
- distinguish odd & even

    x & 1 == 0 //x is even
- Q:Find Repitition:

    Assume an array of 1-max and one more num in [1,max]
    Find the repetition.
    ```cpp
    int FindRepetition(int arr[]) {
        int ans = 0;
        for (int i = 0; i < max+1; i++) {
            ans ^= (arr[i]^(i+1));
        }
        ans ^= (max+1);
        return ans;
    }
    ```
## Tues 03/10/2023
- Q:count "1"s in a integer
    ```cpp
    //让num循环右移与1相与
    int appr1(int num) {
        int ans = 0;
        do {
            ans += num & 1;
            num >>= 1;
        } while (num != 0);
        return ans;
    }
    //让1循环左移与num相与，数表达式true的次数
    int appr2(int num) {
        int ans = 0;
        for (int i = 0; i < 31; ++i) {
            if (num & (1 << i))
                ans++;
        }
        return ans;
    }
    //让num与num-1相与，每次执行ans都++，直到 num & num-1 == 0
    int appr3(int num) {
        int ans = 0;
        while (num != 0) {
            num &= (num - 1);
            ans++;
        }
        return ans;
    }
    ```
- Q:determine whether an integer is n times power of 2 in one single statement
    ```cpp
    //easy one
    bool ans(int num) {
        return !(num & (num-1));
    }
    ```
- Q : print out binary version of a double in (0,1), print 'ERROR' if it can not be presented precisely in 32 numbers
    ```cpp
    //not a hard one, but need to work on strings in cpp
    std::string convert(double num) {
        std::string str = "0.";
        while(num != 0.0) {
            num *= 2.0;
            if(num >= 1.0) {
                str.append("1");
                num -= 1;
            }
            else str.append("0");
        }
        if(str.length() >= 34) {
            str = "ERROR";
        }
        return str;
    }
    ```
- Q : In an array one number shows up for 1 time, others show up k times, find the only number.
- [10/24/23](#102423)
#10/13/2023
- find once number
```cpp
int power(int a, int n) {
    if(n == 0) return 1;
    if(n == 1) return a;
    return a*power(a,n-1);
}

int FindOnceNumber(int vec[], int n, int k) {
    int bits[32];
    for(int i=0; i<32; ++i) {
        bits[i] = 0;
        for(int j=0; j<n; j++){
            bits[i] += ((vec[j] >> i) & 1);
        }
        bits[i] %= k;
        //std::cout << bits[i] << " ";
    }
    int ans = 0;
    for(int i=0; i<32; i++) {
        ans += power(2,i)*bits[i];
    }
    return ans;
}
```
#10/15/2023
- print i~j
```cpp
//ez one
void Printij(int i, int j) {
    std::cout << i << " ";
    if(i == j) return;
    Printij(i+1,j);
}
```
#10/17/2023
- Shell sort
    ```cpp
    void ShellSort(int arr[], int n) {
        //gap为步长
        for(int gap = n >> 1; gap > 0; gap >>= 1) {
            //group为组号，值为0~gap-1
            for(int group = 0; group < gap; group++) {
                //num从group+gap到n，间隔gap，从group+gap开始是为了从第二项开始插入排序，跳过第一项
                for(int num = group+gap; num < n; num += gap) {
                    //插入排序
                    //if (arr[num] < arr[num - gap]) //可以加入这样一个判断，会少运算量
                    int key = arr[num];
                    int j = num-gap;
                    while(key < arr[j] && j >= 0){
                        arr[j+gap] = arr[j];
                        j -= gap;
                    }
                    arr[j+gap] = key;
                }
            }
        }
    }
    ```
- SelectSort
    ```cpp
    void SelectSort(int arr[], int n) {
        for(int i = n-1; i > 0; i--) {
            int max = i;
            for(int j = 0; j < i; j++) {
                if(arr[j] > arr[max]) {
                    max = j;
                }
            }
            swap(arr[max],arr[i]);
        }
    }
    ```
- BubbleSort
    ```cpp
    void BubbleSort(int arr[],int n) {
        for(int i=n-1; i>0; i--) {
            for(int j=0; j<i; j++) {
                if(arr[j] > arr[j+1]) {
                    swap(arr[j],arr[j+1]);
                }
            }
        }
    }
    ```
- Insertsort
    ```cpp
    void InsertSort(int arr[], int n) {
        if(n == 1) return;
        InsertSort(arr,n-1);
        if(arr[n-1] < arr[n-2]) {
            int num = arr[n-1];
            int key = n-2;
            while(num < arr[key] && key >= 0) {
                arr[key+1] = arr[key];
                key--;
            }
            arr[key+1] = num;
        }
    }
    ```
# 10/19/23
- count how many ways to go up n stairs in steps of 1,2 and 3.
    ```cpp
    int stairs(int n) {
        if(n < 0) return 0;
        if(n == 0) return 1;
        return stairs(n-1)+stairs(n-2)+stairs(n-3);
    }
    ```
- search the minimum number of an array 
  input:1,2,3,4,5 2,3,4,5,1 3,4,5,1,2 4,5,1,2,3 5,1,2,3,4
  output: 1
    ```cpp
    int search(int arr[], int n) {
        int begin = 0;
        int end = n-1;
        int mid = begin + (end - begin) >> 1;
        //考虑没有旋转的情况
        if(arr[begin] < arr[end]) return arr[0];
        //循环
        while(begin + 1 < end) {
            mid = begin + ((end - begin) >> 1);
            if(arr[mid] > arr[begin]) {
                begin = mid;
            }
            else {
                end = mid;
            }
        }
        return arr[end];
    }
    ```
- find string in a string array with lots of empty strings
    ```cpp
    int FindString(std::string arr[],int n, std::string str) {
        int begin = 0;
        int end = n-1;
        while(begin < end) {
            int mid = begin + ((end - begin) >> 1);
            int real_mid = mid;
            while(arr[real_mid].empty()) {
                real_mid++;
                if(real_mid == n) {
                    end = begin + ((end - begin) >> 1) - 1;
                    mid = begin + ((end - begin) >> 1);
                    real_mid = mid;
                }
            }
            if(str > arr[real_mid]) {
                begin = real_mid;
            }
            else if(str < arr[real_mid]) {
                end = mid;
            }
            else {
                return real_mid;
            }
        }
        return -1;
    }
    ```
- find the longest ascent subsequence 
    ```cpp
    std::vector<int> FindMaxAscentSeq(std::vector<int> arr) {
        //初始化迭代器变量
        std::pair<std::vector<int>::iterator,std::vector<int>::iterator> ptr(arr.begin(),arr.begin());
        //初始化记录最长子序列迭代�?
        std::pair<std::vector<int>::iterator,std::vector<int>::iterator> max_ptr(arr.begin(),arr.begin());
        //记录当前长度，并初始化为1
        int length = 1;
        //记录最大长度并初始化为0
        int max_length = 0;
        //ptr.second指向end时退出循环
        while(ptr.second != (arr.end()-1)) {
            //当ptr.second小于后一项时让ptr.sencond++,且当前length++
            if(*ptr.second < *(ptr.second + 1)) {
                ptr.second++;
                length++;
            }
            //else，判断
            //如果length > max_length, max_length = length，且记录ptr到max_ptr
            else {
                if(length > max_length) {
                    max_length = length;
                    max_ptr = ptr;
                }
            length = 1;
            //让ptr.first指向之前second的后一项，second指向现在的first
            ptr.first = (ptr.second + 1);
            ptr.second = ptr.first;
            }
            //如果此时second指向end，那么也判断length是否大于max_length，与上述相同
            if(ptr.second == (arr.end()-1)) {
                if(length > max_length) {
                    max_length = length;
                    max_ptr = ptr;
                }
            }
            //返回ptr.first到second之间的vector
        }
        return std::vector<int>(max_ptr.first,max_ptr.second+1);
    }
    ```
- 快速求a的n次幂
    ```cpp
    int power(int a, int n) {
        int ret = a;
        if(n == 0) return 1;
        int ex = 1;
        while((ex<<1) <= n) {
            ret *= ret;
            ex <<= 1;
        }
        return ret*power(a,n-ex);
    }
    ```
#10/21/2023
- 单项扫描分区快速排序
    ```cpp
    /**
    * @brief do the partition part of quicksort in the range [begin,end)
    * @param vec input std::vector<int>
    * @param begin iterator of begin element
    * @param end iterator of next one of the last element
    */
    std::vector<int>::iterator Partition(std::vector<int> &vec, std::vector<int>::iterator begin, std::vector<int>::iterator end) {
        auto pivot = begin;
        auto bg = end - 1;
        auto sp = begin + 1;
        while(sp <= bg) {
            if(*sp <= *pivot) {
                sp++;
            }
            else {
                std::swap(*sp,*bg);
                bg--;
            }
        }
        std::swap(*bg,*pivot);
        return bg;
    }
    /**
     * @brief realization of quicksort algorithm
     * @param vec input std::vector<int>
     * @param begin iterator of begin element
     * @param end iterator of next one of the last element
    */
    void QuickSort(std::vector<int> &vec, std::vector<int>::iterator begin, std::vector<int>::iterator end) {
        if(begin < end-1) {
            auto mid = Partition(vec,begin,end);
            QuickSort(vec,begin,mid);
            QuickSort(vec,mid+1,end);
        }
        else return;
    }
    ```
- 快排的双向扫描分区法
    ```cpp
    std::vector<int>::iterator Partition_2d(std::vector<int> &vec, std::vector<int>::iterator begin, std::vector<int>::iterator end) {
        auto pivot = begin;
        auto left = begin + 1;
        auto right = end - 1;
        while(left <= right) {
            while(left <= right && *left <= *pivot) {
                left++;
            }
            while(left <= right && *right > *pivot) {
                right--;
            }
            if(left < right) {
                std::swap(*left,*right);
            }
        }
        std::swap(*right,*pivot);
        return right;
    }
    ```
- 三指针分区法快排
    ```cpp
    std::pair<std::vector<int>::iterator,std::vector<int>::iterator> Partition_3d(std::vector<int> &vec, std::vector<int>::iterator begin, std::vector<int>::iterator end) {
        auto Next_Less_Pos = begin + 1;
        auto Next_Scan_Pos = Next_Less_Pos;
        auto Next_Bigger_Pos = end - 1;
        auto pivot = begin;
        while(Next_Scan_Pos <= Next_Bigger_Pos) {
            if(*Next_Scan_Pos < *pivot) {
                std::swap(*Next_Less_Pos,*Next_Scan_Pos);
                Next_Less_Pos++;
                Next_Scan_Pos++;
            }
            else if(*Next_Scan_Pos == *pivot) {
                Next_Scan_Pos++;
            }
            else {
                std::swap(*Next_Scan_Pos,*Next_Bigger_Pos);
                Next_Bigger_Pos--;
            }
        }
        Next_Less_Pos--;
        std::swap(*Next_Less_Pos,*pivot);
        return std::pair<std::vector<int>::iterator,std::vector<int>::iterator>(Next_Less_Pos,Next_Bigger_Pos);
    }
    /**
     * @brief realization of quicksort algorithm
     * @param vec input std::vector<int>
     * @param begin iterator of begin element
     * @param end iterator of next one of the last element
    */
    void QuickSort(std::vector<int> &vec, std::vector<int>::iterator begin, std::vector<int>::iterator end) {
        if(begin < end-1) {
            auto mid = Partition_3d(vec,begin,end);
            QuickSort(vec,begin,mid.first);
            QuickSort(vec,mid.second+1,end);
        }
        else return;
    }
    ```
# 10/24/23
- 快排的优点
    - 三点中值法：数组的起始、中间、结尾三个元素作比较得出中间值作为pivot：以双指针分区为例
    ```cpp
    std::vector<int>::iterator Partition_2d(std::vector<int> &vec, std::vector<int>::iterator begin, std::vector<int>::iterator end) {
    //auto pivot = begin;
    auto mid = begin + ((end - begin) >> 1);
    std::vector<int>::iterator pivot;
    if(*begin >= *end && *begin >= *mid) {
        pivot = (*end >= *mid)?end:mid;
    }
    if(*mid >= *end && *mid >= *begin) {
        pivot = (*end >= *begin)?end:begin;
    }
    if(*end >= *begin && *end >= *mid) {
        pivot = (*begin >= *mid)?begin:mid;
    }
    std::swap(*pivot,*begin);
    pivot = begin;
    auto left = begin + 1;
    auto right = end - 1;
    while(left <= right) {
        while(left <= right && *left <= *pivot) {
            left++;
        }
        while(left <= right && *right > *pivot) {
            right--;
        }
        if(left < right) {
            std::swap(*left,*right);
        }
    }
    std::swap(*right,*pivot);
    return right;
    }
    ```
- 绝对中值法
    ```cpp
        //绝对中值法,有点麻烦，写不下去了
    std::vector<int> mid_vector((end-begin)/5+1);
    int gap = 5;
    for(auto it=begin; it<begin+gap; it++) {//5个元素一组
        std::vector<int>::iterator sub_it;
        for(sub_it=it+gap; sub_it < end; sub_it+=gap) {//从每组第二个元素开始插入排序
            if(*sub_it < *(sub_it-gap)) {
                auto key = sub_it;
                auto j = key - gap;
                while(*key < *j && j>=begin) {
                    *(j+gap) = *j;
                    j -= gap;
                }
                *(j+gap) = *key;
            }
        }
    mid_vector.push_back(*(it+((sub_it-it)>>1)));
    std::sort(mid_vector.begin(),mid_vector.end());
    }
    //可以用<algorithm> std::find()
    ```
- 待排序数组元素数小于等于8，用插入排序快于快速排序
- 归并排序
    ```cpp
    void Merge(std::vector<int> &vec,std::vector<int>::iterator begin,std::vector<int>::iterator end) {
        std::vector<int> arr(begin,end);//这里一定记得不要让arr = vec，不然越界！！
        auto left = arr.begin();
        auto new_end = arr.end();
        auto new_mid = left + ((new_end - left) >> 1);
        auto right = new_mid;
        auto ptr = begin;//不要用vec.begin()
        while(left < new_mid && right < new_end) {
            if(*left < *right) {
                *ptr = *left;
                left++;
            }
            else {
                *ptr = *right;
                right++;
            }
            ptr++;
        }
        if(left < new_mid) {
            std::copy(left,new_mid,ptr);
        }
    }

    void MergeSort(std::vector<int> &vec,std::vector<int>::iterator begin,std::vector<int>::iterator end) {
        if(begin + 1 < end) {
            auto mid = begin + ((end - begin) >> 1);
            MergeSort(vec,begin,mid);
            MergeSort(vec,mid,end);
            Merge(vec,begin,end);
        }
    }
    ```
- 把给定整数数组里奇数排在前面，偶数排在后面
    ```cpp
    void Separate_Odd_Even(std::vector<int> &vec) {
        auto sp = vec.begin();
        auto bg = vec.end()-1;
        while(sp <= bg) {
            if((*sp&1)==1){//odd
                sp++;
            }
            else {
                std::swap(*sp,*bg);
                bg--;
            }
        }
    }


    ```
- 以尽量高的效率求出一个乱序数组中按数值顺序的第K个元素
    ```cpp
    std::vector<int>::iterator Partition(std::vector<int> &vec,std::vector<int>::iterator begin,std::vector<int>::iterator end) {
        auto pivot = begin;
        auto sp = begin + 1;
        auto bg = end - 1;
        while(sp <= bg) {
            if(*sp <= *pivot) {
                sp++;
            }
            else {
                std::swap(*sp,*bg);
                bg--;
            }
        }
        std::swap(*pivot,*bg);
        return bg;
    }

    std::vector<int>::iterator SelectK(std::vector<int> &vec,std::vector<int>::iterator begin,std::vector<int>::iterator end,int k) {
        auto mid = Partition(vec,begin,end);
        if(begin + k - 1 == mid) return mid;
        else if(begin + k - 1 > mid) return SelectK(vec,mid+1,end,begin+k-mid-1);
        else return SelectK(vec,begin,mid,k); 
    }

    int Partition(int arr[], int begin, int end) {
    int pivot = begin;
    int sp = begin + 1;
    int bg = end;
    while(sp <= bg) {
        if(arr[sp] <= arr[pivot]) {
            sp++;
        }
        else {
            std::swap(arr[sp],arr[bg]);
            bg--;
        }
    }
    std::swap(arr[bg],arr[pivot]);
    return bg;
    }

    void QuickSort(int arr[], int begin, int end) {
        if(begin < end) {
            int mid = Partition(arr,begin,end);
            QuickSort(arr,begin,mid-1);
            QuickSort(arr,mid+1,end);
        }
    }
    //快速排序后返回第K个值
    int TopK_Sort(int arr[], int length,int k) {
        QuickSort(arr,0,length-1);
        return arr[k-1];
    }

    //冒泡排序K次后输出数组第K-1个值
    int TopK_BubbleSort(int arr[], int length, int k) {
        if(k <= length) {
            for(int i = 0;i < k;i++) {
                for(int j = length - 1;j > i;j--) {
                    if(arr[j] < arr[j-1]) {
                        std::swap(arr[j],arr[j-1]);
                    }
                }
            }
            return arr[k-1];
        }
        return -1;
    }

    void MakeMaxHeap(int arr[], int length) {
        for(int i = length/2 - 1;i >= 0;i--) {
            int j = 2*i + 1;
            int temp = arr[i];
            while(j < length) {
                if(j + 1 < length && arr[j+1] > arr[j]) {
                    j++;
                }
                if(temp < arr[j]) {
                    arr[i] = arr[j];
                    i = j;
                    j = 2*i + 1;
                }
                else break;
            }
            arr[i] = temp;
        }
    }

    int TopK_HeapSort(int arr[], int length, int k) {
        MakeMaxHeap(arr,k);
        for(int i = k;i < length;i++) {
            if(arr[i] < arr[0]) {
                //std::swap(arr[i],arr[0]);
                arr[0] = arr[i];
                int ii = 0;
                int j = 1;
                int temp = arr[ii];
                while(j < k) {
                    if(j + 1 < k && arr[j+1] > arr[j]) {
                        j++;
                    }
                    if(arr[j] > temp) {
                        arr[ii] = arr[j];
                        ii = j;
                        j = 2*ii + 1;
                    }
                    else break;
                }
                arr[ii] = temp;
            }
        }
        return arr[0];
    }

    int TopK_Partition(int arr[], int begin, int end, int k) {
        if(begin <= end) {
            int mid = Partition(arr,begin,end);
            if(mid == k - 1) return arr[mid];
            else if(mid > k - 1) return TopK_Partition(arr,begin,mid-1,k);
            else return TopK_Partition(arr,mid+1,end,k);
        }
        return -1;
    }
    ```
- 数组中出现次数超过半数的元素
    - 可以直接sort后输出arr[length/2],时间复杂度期望O(nlogn)
    - 可以用上例方法输出第length/2个元素，时间复杂度期望O(n)
    - 遍历数组消除
    ```cpp
    int FindMoreThanHalf(std::vector<int> vec) {
    int candidate = vec.front();
    int count = 1;
    for(auto it=vec.begin()+1;it!=vec.end();it++) {
        if(count == 0) {
            candidate = *it;
            count = 1;
            continue;
        }
        if(*it == candidate) {
            count++;
        }
        else {
            count--;
        }
    }
    return candidate;
    }
    ```
- 逆序对个数
    ```cpp
    //使用MergeSort改�?
    int count = 0;
    void Merge(std::vector<int> &vec,std::vector<int>::iterator begin,std::vector<int>::iterator end) {
        std::vector<int> new_vec;
        new_vec.assign(begin,end);
        auto new_begin = new_vec.begin();
        auto new_end = new_vec.end();
        auto new_mid = new_begin + ((new_end-new_begin)>>1);
        auto left = new_begin;
        auto right = new_mid;
        while(left != new_mid && right != new_end) {
            if(*left <= *right) {
                *begin = *left;
                left++;
            }
            else {
                *begin = *right;
                right++;
                count += new_mid - left;
            }
            begin++;
        }
            if(left != new_end) {
                std::copy(left,new_mid,begin);
        }   
    }

    void MergeSort(std::vector<int> &vec,std::vector<int>::iterator begin,std::vector<int>::iterator end) {
        if(begin < end - 1) {
            auto mid = begin + ((end-begin)>>1);
            MergeSort(vec,begin,mid);
            MergeSort(vec,mid,end);
            Merge(vec,begin,end);
        }
    }
    ```
## Thu 12/07/23
- 堆排序
```cpp
void makeMaxHeap(int arr[],int length) {
    for(int i = (length - 1) / 2;i >= 0; i--) {
        int j = 2 * i + 1;
        int temp = arr[i];
        while(j < length) {
            if(j + 1 < length && arr[j + 1] > arr[j]) j++;
            if(temp < arr[j]) {
                arr[i] = arr[j];
                i = j;
                j = 2 * i + 1;
            }
            else break;
        }
        arr[i] = temp;
    }
}
void HeapSort(int arr[], int n) {
    for(int i = n;i > 1;i--) {
        MakeMaxHeap(arr,i);
        std::swap(arr[0],arr[i-1]);
    }
}
```
- 计数排序
```cpp
void CountSort(int arr[], int n) {
    int max = arr[0];
    int min = arr[0];
    for(int i = 0;i < n;i++) {
        if(arr[i] > max) {
            max = arr[i];
        }
        if(arr[i] < min) {
            min = arr[i];
        }
    }
    int range = max - min + 1;
    int new_arr[range] = {0};
    for(int i = 0;i < n;i++) {
        new_arr[arr[i]-min]++;
    }
    int j = 0;
    for(int i = 0;i < range;i++) {
        while(new_arr[i]--) {
            arr[j] = i+min;
            j++;
        }
    }
}
```
- 给定已排序数组arr和k，不重复打印arr中所有相加和为k的不降序二元组
```cpp
    void FindSum10_1(int arr[], int n, int k) {
        for(int i = 0;i < n - 1;i++) {
            for(int j = i + 1;j < n;j++) {
                if(arr[i]+arr[j] == k) {
                    std::cout << arr[i] << " " << arr[j] << std::endl;
                }
            }
        }
    }

    void FindSum10_2(int arr[], int n, int k) {
        for(int i = 0;i < n - 1;i++) {
            int begin = i + 1;
            int end = n - 1;
            while(begin <= end) {
                int mid = begin + ((end - begin) >> 1);
                if(arr[mid] + arr[i]> k) {
                    end = mid - 1;
                }
                else if(arr[mid] + arr[i] < k) {
                    begin = mid + 1;
                }
                else {
                    std::cout << arr[i] << " " << arr[mid] << std::endl;
                    break;
                }
            }
        }
    }

    void FindSum10_3(int arr[], int n, int k) {
        int left = 0;
        int right = n - 1;
        while(left < right) {
            if(arr[left]+arr[right] > k) {
                right--;
            }
            else if(arr[left]+arr[right] < k) {
                left++;
            }
            else {
                std::cout << arr[left] << " " << arr[right] << std::endl;
                left++;
                right--;
            }
        }
    }
```
- 拓展：给定排序数组arr和整数k，不重复打印arr中所有相加和为k的严格升序的三元组
```cpp
void FindSum10_3d(int arr[], int length, int sum) {
    for(int i = 0;i < length - 2;i++) {
        int sumo2 = sum - arr[i];
        int left = i+1;
        int right = length-1;
        while(left < right && arr[i] != arr[left]) {
            if(arr[left] + arr[right] > sumo2) {
                right--;
            }
            else if(arr[left] + arr[right] < sumo2) {
                left++;
            }
            else {
                if(arr[left] != arr[left-1] && arr[right] != arr[right+1] && arr[left] != arr[right]) {
                    std::cout << arr[i] << " " << arr[left] << " " << arr[right] << std::endl;
                }
                left++;
                right--;
            }
        }
    }
}
```
- 给定一个无序数组arr,求出需要排序的最短子数组的长度
```cpp
    int FindMinArrayToSort(int arr[], int length) {
        int max = 0;
        int min = 1e5;
        for(int i = 0;i < length - 1;i++) {
            if(arr[i] > arr[i+1]) {
                if(arr[i] > max) {
                    max = arr[i];
                }
                if(arr[i+1] < min) {
                    min = arr[i+1];
                }
            }
        }
        int begin = -1;
        int end = length - 1;
        for(int i = 0;i < length;i++) {
            if(arr[i] > min && begin == -1) {
                begin = i;
            }
            if(arr[i] > max && end == length - 1) {
                end = i - 1;
            }
        }
        return end - begin + 1;
    }
```
- 把数组排列成最小的数
```cpp
    bool cmp(int a, int b) {
        std::string A;
        std::string B;
        A += std::to_string(a);
        A += std::to_string(b);
        B += std::to_string(b);
        B += std::to_string(a);

        return A<B;
    }

    std::string PrintMinNumber(std::vector<int> vec) {
        if(vec.empty()) return "";
        std::string ret;
        std::sort(vec.begin(),vec.end(),cmp);
        for(int i = 0;i < vec.size();i++) {
            ret += std::to_string(vec[i]);
        }
        return ret;
    }
```
# Sun 12/10/23
- 顺时针打印二维数组
```cpp
    void PrintCircle(std::vector<std::vector<int>> vec,int begin_row, int begin_column, int end_row, int end_column) {
        if(begin_row > end_row || begin_column > end_column) return;
        for(int column = begin_column; column <= end_column;column++) {
            std::cout << vec[begin_row][column] << " ";
        }
        for(int row = begin_row + 1; row <= end_row;row++) {
            std::cout << vec[row][end_column] << " ";
        }
        if(begin_row < end_row) {
            for(int column = end_column - 1; column >= begin_column;column--) {
                std::cout << vec[end_row][column] << " ";
            }
        }
        if (begin_column < end_column) {
            for(int row = end_row - 1; row > begin_row;row--) {
                std::cout << vec[row][begin_column] << " ";
            }
        }
        PrintCircle(vec,begin_row+1,begin_column+1,end_row-1,end_column-1);
    }
```
- 零所在行列清零
```cpp
    void Clear0(int (*arr)[3], int row, int column) {
        int record[row][column] = {0};
        for(int i = 0;i < row;i++) {
            for(int j = 0;j < column;j++) {
                if(arr[i][j] == 0) {
                    record[i][j] = 1;
                }
            }
        }
        for(int i = 0;i < row;i++) {
            for(int j = 0;j < column;j++) {
                if(record[i][j] == 1) {
                    for(int r = 0;r < row;r++) {
                        arr[r][j] = 0;
                    }
                    for(int c = 0;c < column;c++) {
                        arr[i][c] = 0;
                    }
                }
            }
        }
    }
```
- Z型打印二维数组
```cpp
    void printZigZagDiagonal(vector<vector<int>>& matrix) {
        int rows = matrix.size();
        if (rows == 0) {
            cout << "Empty matrix." << endl;
            return;
        }

        int cols = matrix[0].size();

        for (int k = 0; k < rows + cols - 1; ++k) {
            if (k % 2 == 0) { // print from bottom to top
                for (int i = min(k, rows - 1); i >= max(0, k - cols + 1); --i) {//
                    cout << matrix[i][k - i] << " ";
                }
            } else { // print from top to bottom
                for (int i = max(0, k - cols + 1); i <= min(k, rows - 1); ++i) {//0 
                    cout << matrix[i][k - i] << " ";
                }
            }
        }
        cout << endl;
    }
```
## 12/28/23
- 给你两个二进制字符串 a 和 b ，以二进制字符串的形式返回它们的和。
### 代码
#### 方法一
```cpp
class Solution {
public:
    string addBinary(string a, string b) {
        string ans;
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());

        int n = max(a.size(), b.size()), carry = 0;
        for (size_t i = 0; i < n; ++i) {
            carry += i < a.size() ? (a.at(i) == '1') : 0;
            carry += i < b.size() ? (b.at(i) == '1') : 0;
            ans.push_back((carry % 2) ? '1' : '0');
            carry /= 2;
        }

        if (carry) {
            ans.push_back('1');
        }
        reverse(ans.begin(), ans.end());
        return ans;
    }
};
```
#### 方法二
```python
class Solution:
    def addBinary(self, a, b) -> str:
        x, y = int(a, 2), int(b, 2)
        while y:
            answer = x ^ y
            carry = (x & y) << 1
            x, y = answer, carry
        return bin(x)[2:]
```
## 12/29/23
- 给定一个按照升序排列的长度为n的整数数组，以及 q 个查询。

对于每个查询，返回一个元素k的起始位置和终止位置（位置从0开始计数）。

如果数组中不存在该元素，则返回“-1 -1”。
```cpp
/**
 * input:6 3
    1 2 2 3 3 4
    3
    4
    5
   output:3 4
    5 5
    -1 - 1
*/
#include <iostream>

using namespace std;

const int N = 100010;

int main() {
    int n,q;
    int qq[N];
    int k;
    cin >> n >> q;
    for(int i = 0;i < n;i++) {
        cin >> qq[i];
    }
    for(int i = 0; i < q;i++) {
        cin >> k;
        int l = 0;
        int r = n - 1;
        int mid;
        while(l < r) {
            mid = (l + r) >> 1;
            if(qq[mid] < k) l = mid + 1;
            else r = mid;
        }
        if(qq[l] != k) cout << "-1 - 1" << endl;
        else{
            cout << l << " ";
            l = 0;
            r = n - 1;
            int mid;
            while(l < r) {
                mid = (l + r + 1) >> 1;//两种中点
                if(qq[mid] > k) r = mid - 1;
                else l = mid;
            }
            cout << l << endl;
        }
    }
    return 0;
}
```
## 12/30/23
- 高精度运算
```cpp
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
using namespace std;

vector<int> add(vector<int> a, vector<int> b) {
    vector<int> c;
    int t = 0;
    for(int i = 0;i < a.size() || i < b.size();i++) {
        if(i < a.size()) t += a[i];
        if(i < b.size()) t += b[i]; 
        c.push_back(t % 10);
        t /= 10;
    }
    if(t) c.push_back(t);
    return c;
}

// >= return true; < return false
bool cmp(vector<int> a,vector<int> b) {
    if(a.size() != b.size()) return a.size() > b.size();
    for(int i = a.size() - 1 ;i >= 0;i ++) {
        if(a[i] != b[i]) return a[i] > b[i];
    }
    return true;
}

vector<int> minu(vector<int> a, vector<int> b) {
    if(!cmp(a,b)) exit(1);
    vector<int> c;
    int t = 0;
    for(int i = 0;i < a.size();i ++) {
        t += a[i];
        if(i < b.size()) t -= b[i];
        c.push_back((t + 10) % 10);
        if(t < 0) t = -1;
        else t = 0;
    }
    while(c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}

vector<int> mul(vector<int> a, int b) {
    vector<int> c;
    for(int i = 0, t = 0;i < a.size() || t;i ++) {
        if(i < a.size()) t += a[i] * b;
        c.push_back(t % 10);
        t /= 10;
    }
    while(c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}

vector<int> div(vector<int> a, int b) {
    vector<int> c;
    for(int i = a.size() - 1,t = 0;i >= 0;i --) {
        t = 10 * t + a[i];
        c.push_back(t / b);
        t %= b;
    }
    reverse(c.begin(),c.end());
    while(c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}

int main() {
    // string str_a,str_b;
    // cin >> str_a >> str_b;
    // vector<int> a,b;
    // for(int i = str_a.size() - 1;i >= 0;i --) {
    //     a.push_back(str_a[i] - '0');
    // }
    // for(int i = str_b.size() - 1;i >= 0;i --) {
    //     b.push_back(str_b[i] - '0');
    // }
    string str_a;
    int b;
    cin >> str_a >> b;
    vector<int> a;
    for(int i = str_a.size() - 1;i >= 0;i --) {
        a.push_back(str_a[i] - '0');
    }
    vector<int> c = div(a,b);
    for(int i = c.size() - 1;i >= 0;i --) {
        cout << c[i];
    }
    return 0;
}
```
- 输入一个长度为n的整数序列。
接下来再输入m个询问，每个询问输入一对l, r。
对于每个询问，输出原序列中从第l个数到第r个数的和。
```cpp
const int N = 100010;

int main() {
    int n,m;
    int q[N];
    int l,r;
    int sum[N];
    cin >> n >> m;
    for(int i = 0;i < n;i++) {
        cin >> q[i];
        if(i == 0) sum[i] = q[0];
        else sum[i] = sum[i - 1] + q[i];
    }
    while(m--) {
        cin >> l >> r;
        cout << sum[r - 1] - sum[l - 2] << endl;
    }
    return 0;
}
```
- 输入一个n行m列的整数矩阵，再输入q个询问，每个询问包含四个整数x1, y1, x2, y2，表示一个子矩阵的左上角坐标和右下角坐标。

对于每个询问输出子矩阵中所有数的和。
```cpp
const int N = 1010;
int n,m,q;
int s[N][N];
int x1,y1,x2,y2;
int sum[N][N];

int main() {
    cin >> n >> m >> q;
    for(int i = 1;i <= n;i ++) {
        for(int j = 1;j <= m;j ++) {
            cin >> s[i][j];
            sum[i][j] = s[i][j] + sum[i][j - 1] + sum[i - 1][j] - sum[i - 1][j - 1];
        }
    }
    while(q--) {
        cin >> x1 >> y1 >> x2 >> y2;
        cout << sum[x2][y2] - sum[x1 - 1][y2] - sum[x2][y1 - 1] + sum[x1 - 1][y1 - 1] << endl;
    }
    return 0;
}
```
- 输入一个长度为n的整数序列。
接下来输入m个操作，每个操作包含三个整数l, r, c，表示将序列中[l, r]之间的每个数加上c。
请你输出进行完所有操作后的序列。
```cpp
const int N = 100010;
int n,m,l,r,c;
int q[N];
int p[N];

int main() {
    cin >> n >> m;
    for(int i = 1;i <= n;i ++) {
        cin >> q[i];
        p[i] = q[i] - q[i - 1];
    }
    while(m--) {
        cin >> l >> r >> c;
        p[l] += c;
        p[r + 1] -= c;
    }
    for(int i = 1;i <= n;i ++) {
        q[i] = p[i] + q[i - 1];
        cout << q[i] << " ";
    }
    return 0;
}
```
- 输入一个n行m列的整数矩阵，再输入q个操作，每个操作包含五个整数x1, y1, x2, y2, c，其中(x1, y1)和(x2, y2)表示一个子矩阵的左上角坐标和右下角坐标。
每个操作都要将选中的子矩阵中的每个元素的值加上c。
请你将进行完所有操作后的矩阵输出。
```cpp
const int N = 1010;
int n,m,q;
int a[N][N];
int b[N][N];
int x1,y1,x2,y2,c;

void insert(int x1, int y1, int x2, int y2, int c) {
    b[x1][y1] += c;
    b[x2 + 1][y1] -= c;
    b[x1][y2 + 1] -= c;
    b[x2 + 1][y2 + 1] += c;
}

int main() {
    cin >> n >> m >> q;
    for(int i = 1;i <= n;i ++) {
        for(int j = 1;j <= m;j ++) {
            cin >> a[i][j];
            insert(i,j,i,j,a[i][j]);
        }
    }
    while(q--) {
        cin >> x1 >> y1 >> x2 >> y2 >> c;
        insert(x1,y1,x2,y2,c);
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            a[i][j] = b[i][j] + a[i - 1][j] + a[i][j - 1] - a[i - 1][j - 1];
            cout << a[i][j] << " ";
        }
        cout << '\n';
    }
    return 0;
}
```
## 01/02/24
- 最长连续不重复子序列
```cpp
#include <iostream>
using namespace std;

const int N = 1e5 + 10;

int arr[N];

int main() {
    int n;
    cin >> n;
    for(int i = 0;i < n;i ++) {
        cin >> arr[i];
    }
    int max_l = 1;
    int l = 0;
    int r = 1;
    while(r < n) {
        if(arr[r] != arr[r - 1] && r < n) r++;
        else {
            max_l = max(max_l,r - l);
            l = r;
            while(arr[l] == arr[l + 1] && l < n) l++;
            r = l + 1;
        }
    }
    max_l = max(max_l,r - l);
    cout << max_l << endl;
    return 0;
}
```