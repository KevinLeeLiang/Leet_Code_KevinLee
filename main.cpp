#include "include/createTree.h"
#include <iostream>
#include <unordered_map>

using namespace std;

template<typename T>
void print_vector(vector<T> &nums) {
    for (auto num : nums) {
        cout << num << " ";
    }
    cout << endl;
}

template<typename T>
void print_mat(vector<vector<T>> &mat) {
    for (auto row : mat) {
        for (auto p : row) {
            cout << p << " ";
        }
        cout << endl;
    }
}

unordered_map<TreeNode::TreeNode *, int> memo;

int rob(TreeNode::TreeNode *root) {
    if (!root)
        return 0;
    if (memo.count(root)) {
        return memo[root];
    }
    // 访问当前节点
    int robRoot = root->val;
    if (root->left) {
        robRoot += rob(root->left->left) + rob(root->left->right);
    }
    if (root->right) {
        robRoot += rob(root->right->left) + rob(root->right->right);
    }
    // 不访问当前节点
    int notRobRoot = rob(root->left) + rob(root->right);

    int maxMoney = max(robRoot, notRobRoot);
    memo[root] = maxMoney;

    return maxMoney;
}

void rob_test() {
    TreeNode::TreeNode *root = new TreeNode::TreeNode(2);
    root->left = new TreeNode::TreeNode(1);
    root->right = new TreeNode::TreeNode(3);
    root->left->right = new TreeNode::TreeNode(4);
    //  root->right->right = new TreeNode(1);

    cout << "Maximum amount of money: " << rob(root) << endl;
}

#include <vector>

vector<int> countBits(int num) {
    vector<int> result(num + 1, 0);
    for (int i = 1; i <= num; ++i) {
        if (i % 2 == 0) {
            result[i] = result[i / 2];
        } else {
            result[i] = result[i / 2] + 1;
        }
    }
    return result;
}

int countBits_test() {
    int num = 5;
    vector<int> result = countBits(num);

    cout << "Number of 1 bits for numbers from 0 to " << num << ":" << endl;
    for (int i = 0; i <= num; ++i) {
        cout << i << ": " << result[i] << endl;
    }

    return 0;
}

#include "include/reverseVowels.h"

void reverseVowels_test() {
    string str = "leetcode";
    cout << "Original string: " << str << endl;
    reverseVowels_solution res;
    string reversed_str = res.reverseVowels(str);
    cout << "String with reversed vowels: " << reversed_str << endl;
}

#include "include/topKFrequent.h"

void topKFrequent_test() {
    vector<int> nums = {1, 1, 1, 2, 2, 3};
    int k = 2;
    topKFrequent_solution res;
    vector<int> top_k = res.topKFrequent(nums, k);
    cout << "Top " << k << " frequent elements: ";
    for (int num : top_k) {
        cout << num << " ";
    }
    cout << endl;
}

void intersection_test() {
    // 定义两个数组
    std::vector<int> arr1 = {1, 2, 2, 1};
    std::vector<int> arr2 = {2, 2};

    // 定义一个用于存储交集的向量
    std::vector<int> intersection;

    // 首先对两个数组进行排序
    std::sort(arr1.begin(), arr1.end());
    std::sort(arr2.begin(), arr2.end());
    // 使用std::set_intersection函数找到交集
    std::set_intersection(arr1.begin(), arr1.end(), arr2.begin(), arr2.end(),
                          std::back_inserter(intersection));
    unordered_map<int, int> map;
    vector<int> res;
    for (int num : intersection) {
        if (map.find(num) != map.end())
            continue;
        else {
            map[num]++;
            res.push_back(num);
        }
    }
    for (int t : res) {
        std::cout << t << "," << std::endl;
    }
}

#include "include/SummaryRanges.h"

void SummaryRanges_test() {
    SummaryRanges obj;
    obj.addNum(1);
    vector<vector<int>> param_2 = obj.getIntervals(); // [[1,1]]
    for (const auto &interval : param_2) {
        cout << "[" << interval[0] << ", " << interval[1] << "] ";
    }
    cout << endl;

    obj.addNum(3);
    param_2 = obj.getIntervals(); // [[1,1],[3,3]]
    for (const auto &interval : param_2) {
        cout << "[" << interval[0] << ", " << interval[1] << "] ";
    }
    cout << endl;

    obj.addNum(7);
    param_2 = obj.getIntervals(); // [[1,1],[3,3],[7,7]]
    for (const auto &interval : param_2) {
        cout << "[" << interval[0] << ", " << interval[1] << "] ";
    }
    cout << endl;

    obj.addNum(2);
    param_2 = obj.getIntervals(); // [[1,3],[7,7]]
    for (const auto &interval : param_2) {
        cout << "[" << interval[0] << ", " << interval[1] << "] ";
    }
    cout << endl;

    obj.addNum(6);
    param_2 = obj.getIntervals(); // [[1,3],[6,7]]
    for (const auto &interval : param_2) {
        cout << "[" << interval[0] << ", " << interval[1] << "] ";
    }
    cout << endl;

    obj.addNum(8);
    param_2 = obj.getIntervals(); // [[1,3],[6,8]]
    for (const auto &interval : param_2) {
        cout << "[" << interval[0] << ", " << interval[1] << "] ";
    }
    cout << endl;
}

#include "include/maxEnvelopes.h"

void maxEnvelopes_test() {
    vector<vector<int>> envelopes = {{5, 4},
                                     {6, 4},
                                     {6, 7},
                                     {2, 3}};
    // [[5,4],[6,4],[6,7],[2,3]]

    maxEnvelopes::Solution test;
    std::cout << test.maxEnvelopes(envelopes) << std::endl;
}

#include "include/maxSumSubMatrix.h"

void maxSumSubMatrix_test() {
    vector<vector<int>> matrix = {{1, 0,  1},
                                  {0, -2, 3}};
    int k = 2;
    maxSumSubMatrix::Solution test;
    cout << "Maximum sum of submatrix not exceeding " << k << ": "
         << test.maxSumSubmatrix(matrix, k) << endl;
}

#include "include/isPerfectSquare.h"

void isPerfectSquare_test() {
    int num;
    std::cout << "请输入一个数：";
    std::cin >> num;
    isPerfectSquare::solution test;
    if (test.isPerfectSquare(num)) {
        std::cout << num << " 是有效的完全平方数。" << std::endl;
    } else {
        std::cout << num << " 不是有效的完全平方数。" << std::endl;
    }
}

#include "include/largestDivisibleSubset.h"

void largestDivisibleSubset_test() {
    vector<int> nums = {1, 2, 4, 10, 8};
    largestDivisibleSubset::solution test;
    auto res = test.largestDivisibleSubset(nums);
    std::cout << "[";
    for (auto num : res) {
        std::cout << num << " ";
    }
    std::cout << "]" << std::endl;
}

namespace getSum {
    int getSum(int &a, int &b) {
        while (b != 0) {
            unsigned int carry = (unsigned int) (a & b) << 1;
            a = a ^ b;
            b = carry;
        }
        return a;
    }
}; // namespace getSum
void getSum_test() {
    int a = 5, b = 1;
    cout << "Sum of " << a << " and " << b << " is: " << getSum::getSum(a, b)
         << endl;
}

namespace kSmallestPairs {
    vector<vector<int>> kSmallestPairs(vector<int> &nums1, vector<int> &nums2,
                                       int k) {
        auto cmp = [&nums1, &nums2](const pair<int, int> &a,
                                    const pair<int, int> &b) {
            return nums1[a.first] + nums2[a.second] > nums1[b.first] + nums2[b.second];
        };
        priority_queue<pair<int, int>, vector<pair<int, int >>, decltype(cmp)> pq(cmp);
        int m = nums1.size();
        int n = nums2.size();
        vector<vector<int>> ans;
        for (int i = 0; i < min(k, m); ++i) {
            pq.emplace(i, 0);
        }
        while (k-- > 0 && !pq.empty()) {
            auto[x, y] = pq.top();
            pq.pop();
            ans.emplace_back(initializer_list<int>{nums1[x], nums2[y]});
            if (y + 1 < n) {
                pq.emplace(x, y + 1);
            }
        }
        return ans;
    }
} // namespace kSmallestPairs

void kSmallestPairs_test() {
    vector<int> nums1 = {1, 2, 4, 5, 6};
    vector<int> nums2 = {3, 5, 7, 9};
    int k = 3;
    vector<vector<int>> pairs = kSmallestPairs::kSmallestPairs(nums1, nums2, k);

    cout << "K smallest pairs:" << endl;
    for (const auto &pair : pairs) {
        cout << pair[0] << " " << pair[1] << endl;
    }
}

namespace wiggleMaxLength {
    int wiggleMaxLength(vector<int> &nums) {
        int n = nums.size();
        if (n < 2) {
            return n;
        }
        vector<int> up(n), down(n);
        up[0] = down[0] = 1;
        for (int i = 1; i < n; ++i) {
            if (nums[i] > nums[i - 1]) {
                up[i] = max(up[i - 1], down[i - 1] + 1);
                down[i] = down[i - 1];
            } else if (nums[i] < nums[i - 1]) {
                up[i] = up[i - 1];
                down[i] = max(down[i - 1], up[i - 1] + 1);
            } else {
                up[i] = up[i - 1];
                down[i] = down[i - 1];
            }
        }
        return max(up[n - 1], down[n - 1]);
    }
} // namespace wiggleMaxLength

void wiggleMaxLength_test() {
    // vector<int>nums = {1,17,5,10,13,15,10,5,16,8};
    vector<int> nums = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    cout << wiggleMaxLength::wiggleMaxLength(nums) << endl;
}

namespace canConstruct {
    bool canConstruct(string ransomNote, string magazine) {
        unordered_map<char, int> map;
        for (char c : magazine) {
            map[c]++;
        }
        for (char c : ransomNote) {
            if (map[c] > 0) {
                map[c]--;
            } else {
                return false;
            }
        }
        return true;
    }
} // namespace canConstruct

void canConstruct_test() {
    string ransomNote = "aabbcc";
    string magazine = "abcabcabc";

    if (canConstruct::canConstruct(ransomNote, magazine)) {
        cout << "可以构成赎金信" << endl;
    } else {
        cout << "无法构成赎金信" << endl;
    }
}

#include "include/ShuffleArray.h"

void shuffle_test() {
    std::vector<int> nums = {1, 2, 3};
    Shuffle::ShuffleArray obj(nums);
    std::vector<int> param_1 = obj.reset();
    std::vector<int> param_2 = obj.shuffle();
    for (int num : param_1) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
    for (int num : param_2) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

#include "include/NestedInteger.h"

void NestedInteger_test() {
    NestedInteger::Solution sol;
    string input = "[123,[456,[789]]]";
    NestedInteger::NestedInteger result = sol.deserialize(input);
    std::cout << result.getInteger() << std::endl;
    // Output the result or further process it
}

namespace lexicalOrder {

    void lexicalOrderHelper(int current, int n, vector<int> &result) {
        if (current > n) {
            return;
        }

        result.push_back(current);

        // 在当前数字后面添加一个数字，生成下一个数字
        for (int i = 0; i <= 9; ++i) {
            if (10 * current + i > n) {
                return;
            }
            lexicalOrderHelper(10 * current + i, n, result);
        }
    }

    vector<int> lexicalOrder(int n) {
        vector<int> result;
        for (int i = 1; i <= 9; ++i) {
            lexicalOrderHelper(i, n, result);
        }
        return result;
    }
}; // namespace lexicalOrder

void lexicalOrder_test() {
    int n = 13;
    vector<int> result = lexicalOrder::lexicalOrder(n);

    cout << "Lexical order from 1 to " << n << ":" << endl;
    for (int num : result) {
        cout << num << " ";
    }
    cout << endl;
}

namespace firstUniqChar {
    int firstUniqChar(const string &s) {
        unordered_map<char, int> map;
        for (char c : s) {
            map[c]++;
        }
        // 第二次遍历，找到第一个出现次数为 1 的字符
        for (int i = 0; i < s.length(); ++i) {
            if (map[s[i]] == 1) {
                return i;
            }
        }

        return -1; // 如果没有唯一字符，返回 -1
    }
} // namespace firstUniqChar

void firstUniqChar_test() {
    string s = "leetcode";
    int index = firstUniqChar::firstUniqChar(s);

    if (index != -1) {
        cout << "The first unique character is '" << s[index] << "' at index "
             << index << "." << endl;
    } else {
        cout << "There is no unique character in the string." << endl;
    }
}

namespace lastRemaining {
    int lastRemaining(int n) {
        int a1 = 1;
        int k = 0, cnt = n, step = 1;
        while (cnt > 1) {
            if (k % 2 == 0) { // 正向
                a1 = a1 + step;
            } else { // 反向
                a1 = (cnt % 2 == 0) ? a1 : a1 + step;
            }
            k++;
            cnt = cnt >> 1;
            step = step << 1;
        }
        return a1;
    }
} // namespace lastRemaining
void lastRemaining_test() {
    int n = 9;
    cout << lastRemaining::lastRemaining(n) << endl;
}

namespace isRectangleCover {
    typedef pair<int, int> Point;

    bool isRectangleCover(vector<vector<int>> &rectangles) {
        long area = 0;
        int minX = rectangles[0][0], minY = rectangles[0][1], maxX = rectangles[0][2],
                maxY = rectangles[0][3];
        map<Point, int> cnt;
        for (auto &rect : rectangles) {
            int x = rect[0], y = rect[1], a = rect[2], b = rect[3];
            area += (long) (a - x) * (b - y);

            minX = min(minX, x);
            minY = min(minY, y);
            maxX = max(maxX, a);
            maxY = max(maxY, b);

            Point point1({x, y});
            Point point2({x, b});
            Point point3({a, y});
            Point point4({a, b});

            cnt[point1] += 1;
            cnt[point2] += 1;
            cnt[point3] += 1;
            cnt[point4] += 1;
        }

        Point pointMinMin({minX, minY});
        Point pointMinMax({minX, maxY});
        Point pointMaxMin({maxX, minY});
        Point pointMaxMax({maxX, maxY});
        if (area != (long long) (maxX - minX) * (maxY - minY) ||
            !cnt.count(pointMinMin) || !cnt.count(pointMinMax) ||
            !cnt.count(pointMaxMin) || !cnt.count(pointMaxMax)) {
            return false;
        }

        cnt.erase(pointMinMin);
        cnt.erase(pointMinMax);
        cnt.erase(pointMaxMin);
        cnt.erase(pointMaxMax);

        for (auto &entry : cnt) {
            int value = entry.second;
            if (value != 2 && value != 4) {
                return false;
            }
        }
        return true;
    }
} // namespace isRectangleCover

void isRectangleCover_test() {
    vector<vector<int>> rectangles = {
            {1, 1, 3, 3},
            {3, 1, 4, 2},
            {3, 2, 4, 4},
            {1, 3, 2, 4},
            {2, 3, 3, 4}};
    if (isRectangleCover::isRectangleCover(rectangles)) {
        cout << "能够组成完美矩形" << endl;
    } else {
        cout << "不能组成完美矩形" << endl;
    }
    rectangles = {{1, 1, 2, 3},
                  {1, 3, 2, 4},
                  {3, 1, 4, 2},
                  {3, 2, 4, 4}};
    if (isRectangleCover::isRectangleCover(rectangles)) {
        cout << "能够组成完美矩形" << endl;
    } else {
        cout << "不能组成完美矩形" << endl;
    }
    rectangles = {{0, 0, 1, 1},
                  {0, 1, 3, 2},
                  {1, 0, 2, 2}};
    if (isRectangleCover::isRectangleCover(rectangles)) {
        cout << "能够组成完美矩形" << endl;
    } else {
        cout << "不能组成完美矩形" << endl;
    }
}

#include <stack>

namespace decodeString {
    string getDigits(string &s, size_t &ptr) {
        string ret = "";
        while (isdigit(s[ptr])) {
            ret.push_back(s[ptr++]);
        }
        return ret;
    }

    string getString(vector<string> &v) {
        string ret;
        for (const auto &s : v) {
            ret += s;
        }
        return ret;
    }

    string decodeString(string s) {
        vector<string> stk;
        size_t ptr = 0;
        while (ptr < s.size()) {
            char cur = s[ptr];
            if (isdigit(cur)) {
                string digits = getDigits(s, ptr);
                stk.push_back(digits);
            } else if (isalpha(cur) || cur == '[') {
                stk.push_back(string(1, s[ptr]));
                ptr++;
            } else {
                ptr++;
                vector<string> sub;
                while (stk.back() != "[") {
                    sub.push_back(stk.back());
                    stk.pop_back();
                }
                reverse(sub.begin(), sub.end());
                stk.pop_back();
                // 此时栈顶为当前 sub 对应的字符串应该出现的次数
                int repTime = stoi(stk.back());
                stk.pop_back();
                string t, o = getString(sub);
                while (repTime--)
                    t += o;
                stk.push_back(t);
            }
        }
        return getString(stk);
    }
} // namespace decodeString
void decodeString_test() {
    string s = "3[a2[c]]";
    cout << decodeString::decodeString(s) << endl;
}

namespace longestSubstring {
    int dfs(const string &s, int l, int r, int k) {
        vector<int> cnt(26, 0);
        for (int i = l; i <= r; i++) {
            cnt[s[i] - 'a']++;
        }

        char split = 0;
        for (int i = 0; i < 26; i++) {
            if (cnt[i] > 0 && cnt[i] < k) {
                split = i + 'a';
                break;
            }
        }
        if (split == 0) {
            return r - l + 1;
        }

        int i = l;
        int ret = 0;
        while (i <= r) {
            while (i <= r && s[i] == split) {
                i++;
            }
            if (i > r) {
                break;
            }
            int start = i;
            while (i <= r && s[i] != split) {
                i++;
            }

            int length = dfs(s, start, i - 1, k);
            ret = max(ret, length);
        }
        return ret;
    }

    int longestSubstring(string s, int k) {
        int n = s.length();
        return dfs(s, 0, n - 1, k);
    }
} // namespace longestSubstring

void longestSubstring_test() {
    string s = "aaabb";
    int k = 3;
    cout << "字符串 " << s << " 子串中的每一字符出现次数都不少于 " << k
         << " 这一子串的长度 " << longestSubstring::longestSubstring(s, k)
         << endl;
    s = "ababbc";
    k = 2;
    cout << "字符串 " << s << " 子串中的每一字符出现次数都不少于 " << k
         << " 这一子串的长度 " << longestSubstring::longestSubstring(s, k)
         << endl;
}

namespace maxRotateFunction {
    int maxRotateFunction(vector<int> &nums) {
        int f = 0, n = nums.size();
        int numSum = accumulate(nums.begin(), nums.end(), 0);
        for (int i = 0; i < n; i++) {
            f += i * nums[i];
        }
        int res = f;
        for (int i = n - 1; i > 0; i--) {
            f += numSum - n * nums[i];
            res = max(res, f);
        }
        return res;
    }
} // namespace maxRotateFunction

void maxRotateFunction_test() {
    cout << atan(2.908 / 5.92) * 180 / M_PI << endl;
    vector<int> nums = {4, 3, 2, 6};
    cout << "nums 输入最大值：" << maxRotateFunction::maxRotateFunction(nums)
         << endl;
}

namespace findNthDigit {
    int findNthDigit(int n) {
        int d = 1, count = 9;
        while (n > (long) d * count) {
            n -= d * count;
            d++;
            count *= 10;
        }
        int index = n - 1;
        int start = (int) pow(10, d - 1);
        int num = start + index / d;
        int digitIndex = index % d;
        int digit = (num / (int) (pow(10, d - digitIndex - 1))) % 10;
        return digit;
    }
} // namespace findNthDigit

void findNthDigit_test() {
    int n = 15;
    cout << "第" << n << "位数字为：" << findNthDigit::findNthDigit(n) << endl;
}

namespace removeKdigits {
    string removeKdigits(string num, int k) {
        vector<char> stk;
        for (auto &digit : num) {
            while (stk.size() > 0 && stk.back() > digit && k) {
                stk.pop_back();
                k -= 1;
            }
            stk.push_back(digit);
        }

        for (; k > 0; --k) {
            stk.pop_back();
        }

        string ans = "";
        bool isLeadingZero = true;
        for (auto &digit : stk) {
            if (isLeadingZero && digit == '0') {
                continue;
            }
            isLeadingZero = false;
            ans += digit;
        }
        return ans == "" ? "0" : ans;
    }
} // namespace removeKdigits

void removeKdigits_test() {
    string num = "1432219";
    int k = 3;
    std::cout << "num:" << num << "remove " << k << " 位数字后的最小数字是 "
              << removeKdigits::removeKdigits(num, k) << endl;
}

namespace canCross {
    bool canCross(vector<int> &stones) {
        int n = stones.size();
        vector<vector<int>> dp(n, vector<int>(n));
        dp[0][0] = true;
        for (int i = 1; i < n; ++i) {
            if (stones[i] - stones[i - 1] > i) {
                return false;
            }
        }
        for (int i = 1; i < n; i++) {
            for (int j = i - 1; j >= 0; --j) {
                int k = stones[i] - stones[j];
                if (k > i + 1) {
                    break;
                }
                dp[i][k] = dp[j][k - 1] || dp[j][k] || dp[j][k + 1];
                if (i == n - 1 && dp[i][k])
                    return true;
            }
        }
        return false;
    }
}; // namespace canCross

void canCross_test() {
    std::vector<int> stones = {0, 1, 3, 5, 6, 8, 12, 17}; // 石头的位置
    bool canCrossRiver = canCross::canCross(stones);

    if (canCrossRiver) {
        std::cout << "青蛙可以成功过河！" << std::endl;
    } else {
        std::cout << "青蛙无法成功过河。" << std::endl;
    }
}

namespace fizzBuzz {
    vector<string> fizzBuzz(int n) {
        vector<string> res;
        for (int i = 1; i <= n; ++i) {
            res.push_back((i % 3 == 0 ? (i % 5 == 0 ? "FizzBuzz" : "Fizz")
                                      : (i % 5 == 0 ? "Buzz" : to_string(i))));
        }
        return res;
    }
} // namespace fizzBuzz

void fizzBuzz_test() {
    int n = 15;
    vector<string> res = fizzBuzz::fizzBuzz(n);
    for (int i = 0; i < res.size(); ++i) {
        if (i < res.size() - 1)
            std::cout << res[i] << ", ";
        else
            std::cout << res[i] << std::endl;
    }
}

namespace numberOfArithmeticSlices {
    int numberOfArithmeticSlices(vector<int> &nums) {
        int n = nums.size();
        if (n < 3)
            return 0;
        vector<int> dp(n, 0);

        for (int i = 2; i < n; ++i) {
            if (nums[i - 2] - nums[i - 1] == nums[i - 1] - nums[i])
                dp[i] = dp[i - 1] + 1;
            //      else
            //        dp[i] = dp[i - 1];
        }
        int result = 0;
        for (int i = 0; i < n; ++i) {
            result += dp[i];
        }

        return result;
    }
} // namespace numberOfArithmeticSlices

void numberOfArithmeticSlices_test() {
    vector<int> nums = {1, 2, 3, 8, 9, 10};
    nums = {1, 2, 3, 4};
    std::cout << numberOfArithmeticSlices::numberOfArithmeticSlices(nums)
              << std::endl;
}

namespace thirdMax {
    int thirdMax(vector<int> &nums) {
        sort(nums.begin(), nums.end(), greater<>());
        for (int i = 1, diff = 1; i < nums.size(); ++i) {
            if (nums[i] != nums[i - 1] && ++diff == 3) { // 此时 nums[i] 就是第三大的数
                return nums[i];
            }
        }
        return nums[0];
    }
} // namespace thirdMax

void thirdMax_test() {
    vector<int> nums = {1, 2, 2};
    std::cout << thirdMax::thirdMax(nums) << std::endl;
}

namespace addString {
    string addStrings(string num1, string num2) {
        int i = num1.length() - 1, j = num2.length() - 1, add = 0;
        string ans = "";
        while (i >= 0 || j >= 0 || add != 0) {
            int x = i >= 0 ? num1[i] - '0' : 0;
            int y = j >= 0 ? num2[j] - '0' : 0;
            int result = x + y + add;
            ans.push_back('0' + result % 10);
            add = result / 10;
            i -= 1;
            j -= 1;
        }
        // 计算完以后的答案需要翻转过来
        reverse(ans.begin(), ans.end());
        return ans;
    }
} // namespace addString
void addStrings_test() {
    string num1 = "456";
    string num2 = "77";
    std::cout << num1 << " add " << num2 << " is "
              << addString::addStrings(num1, num2) << std::endl;
}

namespace canPartition {
    bool canPartition(vector<int> &nums) {
        int sum = std::accumulate(nums.begin(), nums.end(), 0);
        if (sum % 2 != 0) {
            return false; // 如果数组元素和为奇数，则无法分割成两个和相等的子集
        }

        int target = sum / 2;
        std::vector<bool> dp(target + 1, false);
        dp[0] = true;

        for (int num : nums) {
            for (int i = target; i >= num; --i) {
                dp[i] = dp[i] || dp[i - num];
            }
        }

        return dp[target];
    }
} // namespace canPartition

void canPartition_test() {
    std::vector<int> nums = {1, 5, 11, 5};
    bool canBePartitioned = canPartition::canPartition(nums);

    if (canBePartitioned) {
        std::cout << "数组可以被分割成两个和相等的子集。" << std::endl;
    } else {
        std::cout << "数组无法被分割成两个和相等的子集。" << std::endl;
    }
}

namespace countBattleships {
    int countBattleships(vector<vector<char>> &board) {
        int ans = 0;
        int rows = board.size();
        int cols = board[0].size();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (board[i][j] == 'X') {
                    if ((i == 0 || board[i - 1][j] != 'X') &&
                        (j == 0 || board[i][j - 1] != 'X')) {
                        ans++;
                    }
                }
            }
        }
        return ans;
    }
} // namespace countBattleships

void countBattleships_test() {
    std::vector<std::vector<char>> board = {
            {'X', '.', '.', 'X'},
            {'.', '.', '.', 'X'},
            {'.', '.', '.', 'X'}};

    std::cout << "Number of battleships: "
              << countBattleships::countBattleships(board) << std::endl;
}

namespace minMutation {
    int minMutation(string startGene, string endGene, vector<string> &bank) {
        unordered_set<string> bank_set(bank.begin(), bank.end());
        if (!bank_set.count(endGene))
            return -1;
        queue<pair<string, int>> q;
        q.push({startGene, 0});
        char genes[] = {'A', 'C', 'G', 'T'};
        while (!q.empty()) {
            auto[current, steps] = q.front();
            q.pop();
            if (current == endGene)
                return steps;
            for (int i = 0; i < current.size(); ++i) {
                char originalChar = current[i];
                for (char gene : genes) {
                    if (gene != originalChar) {
                        current[i] = gene;
                        if (bank_set.count(current)) {
                            q.push({current, steps + 1});
                            bank_set.erase(current);
                        }
                    }
                }
                current[i] = originalChar;
            }
        }
        return -1; // 无法达到目标基因序列
    }
} // namespace minMutation

void minMutation_test() {
    string start = "AACCGGTT";
    string end = "AAACGGTA";
    vector<string> bank = {"AACCGGTA", "AACCGCTA", "AAACGGTA"};

    int result = minMutation::minMutation(start, end, bank);
    if (result != -1) {
        cout << "最小基因变化次数: " << result << endl;
    } else {
        cout << "无法达到目标基因序列" << endl;
    }
}

namespace countSegments {
    int countSegments(string s) {
        int ans = 0;
        string tmp = "";
        for (size_t i = 0; i < s.size(); ++i) {
            char c = s[i];
            if (c == ' ') {
                if (tmp != "") {
                    ans++;
                    tmp = "";
                } else {
                    continue;
                }
            } else {
                if (i == s.size() - 1)
                    ans++;
                else
                    tmp += c;
            }
        }
        return ans;
    }
} // namespace countSegments

void countSegments_test() {
    string s = "  Hello, my name is John";
    cout << s << " 共有：" << countSegments::countSegments(s) << endl;
}

namespace levelOrder {
    class Node {
    public:
        int val;
        vector<Node *> children;

        Node() {}

        Node(int _val) { val = _val; }

        Node(int _val, vector<Node *> _children) {
            val = _val;
            children = _children;
        }
    };

    vector<vector<int>> levelOrder(Node *root) {
        if (root == nullptr)
            return {};
        vector<vector<int>> ans;
        queue<Node *> q;
        q.push(root);

        while (!q.empty()) {
            int size = q.size();
            vector<int> tmp;
            for (int i = 0; i < size; ++i) {
                auto p = q.front();
                tmp.push_back(p->val);
                for (auto t : p->children) {
                    q.push(t);
                }
                q.pop();
            }
            ans.push_back(tmp);
        }
        return ans;
    }
} // namespace levelOrder
void levelOrder_test() {
    levelOrder::Node *root = new levelOrder::Node(1);
    root->children.push_back(new levelOrder::Node(2));
    root->children.push_back(new levelOrder::Node(3));
    root->children.push_back(new levelOrder::Node(4));
    root->children.push_back(new levelOrder::Node(5));
    root->children[1]->children.push_back(new levelOrder::Node(6));
    root->children[1]->children.push_back(new levelOrder::Node(7));
    root->children[1]->children[1]->children.push_back(new levelOrder::Node(11));
    root->children[1]->children[1]->children[0]->children.push_back(
            new levelOrder::Node(14));
    root->children[2]->children.push_back(new levelOrder::Node(8));
    root->children[2]->children[0]->children.push_back(new levelOrder::Node(12));
    root->children[3]->children.push_back(new levelOrder::Node(9));
    root->children[3]->children.push_back(new levelOrder::Node(10));
    root->children[3]->children[0]->children.push_back(new levelOrder::Node(13));
    auto ans = levelOrder::levelOrder(root);
    for (auto list : ans) {
        for (auto p : list) {
            cout << p << " ";
        }
        cout << endl;
    }
}

namespace eraseOverlapIntervals {
    int eraseOverlapIntervals(vector<vector<int>> &intervals) {
        if (intervals.empty()) {
            return 0;
        }

        sort(intervals.begin(), intervals.end(),
             [](const auto &u, const auto &v) { return u[1] < v[1]; });
        int n = intervals.size();
        int right = intervals[0][1];
        int ans = 1;
        for (int i = 1; i < n; ++i) {
            if (intervals[i][0] >= right) {
                ++ans;
                right = intervals[i][1];
            }
        }
        return n - ans;
    }
} // namespace eraseOverlapIntervals

void eraseOverlapIntervals_test() {
    vector<vector<int>> intervals = {{1, 2},
                                     {2, 3},
                                     {3, 4},
                                     {1, 3}};
    cout << "移除 " << eraseOverlapIntervals::eraseOverlapIntervals(intervals)
         << " 来使剩下的区间没有重叠" << endl;
}

namespace findRightInterval {
    vector<int> findRightInterval(vector<vector<int>> &intervals) {
        vector<pair<int, int>> startIntervals;
        int n = intervals.size();
        for (int i = 0; i < n; i++) {
            startIntervals.emplace_back(intervals[i][0], i);
        }
        sort(startIntervals.begin(), startIntervals.end());

        vector<int> ans(n, -1);
        for (int i = 0; i < n; i++) {
            auto it = lower_bound(startIntervals.begin(), startIntervals.end(),
                                  make_pair(intervals[i][1], 0));
            if (it != startIntervals.end()) {
                ans[i] = it->second;
            }
        }
        return ans;
    }
} // namespace findRightInterval

void findRightInterval_test() {
    vector<vector<int>> intervals = {{3, 4},
                                     {2, 3},
                                     {1, 2}};
    auto res = findRightInterval::findRightInterval(intervals);
    for (auto tmp : res) {
        cout << tmp << " ";
    }
    cout << endl;
}

TreeNode::TreeNode *create_treenode(vector<int> tree_vals) {
    auto *tree = TreeNode::createTree(tree_vals);
    //	cout << tree->val << endl;
    return tree;
}

TreeNode::TreeNode *create_treenode(vector<int> tree_vals, bool is_include_zero) {
    if (is_include_zero) {
        auto *tree = TreeNode::createTree2(tree_vals);
        return tree;
    } else {
        auto *tree = TreeNode::createTree(tree_vals);
        //	cout << tree->val << endl;
        return tree;
    }
}

namespace pathSum {
    unordered_map<long long, int> prefix;

    int dfs(TreeNode::TreeNode *root, long long curr, int targetSum) {
        if (!root) {
            return 0;
        }

        int ret = 0;
        curr += root->val;
        if (prefix.count(curr - targetSum)) {
            ret = prefix[curr - targetSum];
        }

        prefix[curr]++;
        ret += dfs(root->left, curr, targetSum);
        ret += dfs(root->right, curr, targetSum);
        prefix[curr]--;

        return ret;
    }

    int pathSum(TreeNode::TreeNode *root, int targetSum) {
        prefix[0] = 1;
        return dfs(root, 0, targetSum);
    }
} // namespace pathSum

void pathSum_test() {
    vector<int> tree_vals = {10, 5, -3, 3, 2, 0, 11, 3, -2, 0, 1};
    auto *tree = create_treenode(tree_vals);
    int targetSum = 8;
    cout << pathSum::pathSum(tree, targetSum) << endl;
}

namespace findKthNumber {
    int calculateSteps(int n, long long curr, long long next) {
        int steps = 0;
        while (curr <= n) {
            steps += min((long long) n + 1, next) - curr;
            curr *= 10;
            next *= 10;
        }
        return steps;
    }

    int findKthNumber(int n, int k) {
        int curr = 1;
        k--; // 因为我们是从1开始的，所以先减去1

        while (k > 0) {
            int steps = calculateSteps(n, curr, curr + 1);
            if (steps <= k) {
                // 如果当前前缀下的数字数量小于等于k，跳到下一个前缀
                curr += 1;
                k -= steps;
            } else {
                // 如果当前前缀下的数字数量大于k，深入到该前缀的下一层
                curr *= 10;
                k -= 1;
            }
        }

        return curr;
    }
} // namespace findKthNumber
void findKthNumber_test() {
    int n = 13;
    int k = 3;

    int result = findKthNumber::findKthNumber(n, k);
    cout << "字典序的第" << k << "小数字是: " << result << endl;
}

namespace arrangeCoins {
    int arrangeCoins(int n) {
        int left = 1, right = n;
        while (left < right) {
            int mid = (right - left + 1) / 2 + left;
            if ((long long) mid * (mid + 1) <= (long long) 2 * n) {
                left = mid;
            } else {
                right = mid - 1;
            }
        }
        return left;
    }
} // namespace arrangeCoins

void arrangeCoins_test() {
    int n = 5;
    cout << "给你一个数字 " << n << " ，计算并返回可形成完整阶梯行的总行数为："
         << arrangeCoins::arrangeCoins(n) << endl;
}

namespace findDuplicates {
    vector<int> findDuplicates(vector<int> &nums) {
        int n = nums.size();
        vector<int> ans;
        for (int i = 0; i < n; ++i) {
            int x = abs(nums[i]);
            if (nums[x - 1] > 0) {
                nums[x - 1] = -nums[x - 1];
            } else {
                ans.push_back(x);
            }
        }
        return ans;
    }
} // namespace findDuplicates
void findDuplicates_test() {
    vector<int> nums = {4, 3, 2, 7, 8, 2, 3, 1};
    vector<int> ans = findDuplicates::findDuplicates(nums);
    cout << "数组中重复的数据 [";
    for (auto p : ans) {
        cout << " " << p;
    }
    cout << " ]" << endl;
}

namespace compress {
    int compress(vector<char> &chars) {
        if (chars.size() == 0)
            return 0;
        //    if (chars.size() == 1)
        //      return 1;
        char tmp = chars[0];
        string ans;
        int count = 1;
        for (int i = 1; i < chars.size(); ++i) {
            if (tmp == chars[i]) {
                count++;
                if (i == chars.size() - 1) {
                    ans += tmp;
                    if (count > 1)
                        ans += to_string(count);
                    count = 0;
                }
            } else {
                ans += tmp;
                if (count > 1)
                    ans += to_string(count);
                count = 1;
                tmp = chars[i];
            }
        }
        if (count != 0) {
            ans += tmp;
            if (count > 1)
                ans += to_string(count);
        }
        chars.assign(ans.begin(), ans.end());
        return ans.size();
    }
} // namespace compress

void compress_test() {
    vector<char> chars = {'a', 'a', 'b', 'b', 'c', 'c', 'c'};
    string str1(chars.begin(), chars.end());
    cout << "chars: " << string(chars.begin(), chars.end())
         << " 压缩后的数组的新长度：" << compress::compress(chars)
         << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
    chars = {'a'};
    cout << "chars: " << string(chars.begin(), chars.end())
         << " 压缩后的数组的新长度：" << compress::compress(chars)
         << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
    chars = {'a', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b'};
    cout << "chars: " << string(chars.begin(), chars.end())
         << " 压缩后的数组的新长度：" << compress::compress(chars)
         << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
    chars = {'a', 'a', 'a', 'b', 'b', 'a', 'a'};
    cout << "chars: " << string(chars.begin(), chars.end())
         << " 压缩后的数组的新长度：" << compress::compress(chars)
         << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
    chars = {'a', 'b', 'c'};
    cout << "chars: " << string(chars.begin(), chars.end())
         << " 压缩后的数组的新长度：" << compress::compress(chars)
         << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
}

namespace numberOfBoomerangs {
    int numberOfBoomerangs(vector<vector<int>> &points) {
        int ans = 0;
        for (auto &p : points) {
            unordered_map<int, int> cnt;
            for (auto &q : points) {
                int dis = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]);
                ++cnt[dis];
            }
            for (auto &[_, m] : cnt) {
                ans += m * (m - 1);
            }
        }
        return ans;
    }
}; // namespace numberOfBoomerangs

void numberOfBoomerangs_test() {
    vector<vector<int>> points = {{0, 0},
                                  {1, 0},
                                  {2, 0}};
    cout << "points number of boomerangs is "
         << numberOfBoomerangs::numberOfBoomerangs(points) << endl;
}

namespace findDisappearedNumbers {
    vector<int> findDisappearedNumbers(vector<int> &nums) {
        int n = nums.size();
        for (auto &num : nums) {
            int x = (num - 1) % n;
            nums[x] += n;
        }
        vector<int> ret;
        for (int i = 0; i < n; i++) {
            if (nums[i] <= n) {
                ret.push_back(i + 1);
            }
        }
        return ret;
    }
} // namespace findDisappearedNumbers
void findDisappearedNumbers_test() {
    vector<int> nums = {4, 3, 2, 7, 8, 2, 3, 1};
    auto ans = findDisappearedNumbers::findDisappearedNumbers(nums);
    cout << "消失的字数字：";
    for (auto p : ans) {
        cout << " " << p;
    }
    cout << endl;
}

namespace SerializingAndDeserializingForBinaryTrees {
    vector<string> split(const string &str, char dec) {
        int pos = 0;
        int start = 0;
        vector<string> res;
        while (pos < str.size()) {
            while (pos < str.size() && str[pos] == dec) {
                pos++;
            }
            start = pos;
            while (pos < str.size() && str[pos] != dec) {
                pos++;
            }
            if (start < str.size()) {
                res.emplace_back(str.substr(start, pos - start));
            }
        }
        return res;
    }

    void postOrder(TreeNode::TreeNode *root, vector<int> &arr) {
        if (root == nullptr) {
            return;
        }
        postOrder(root->left, arr);
        postOrder(root->right, arr);
        arr.emplace_back(root->val);
    }

    TreeNode::TreeNode *construct(int lower, int upper, stack<int> &st) {
        if (st.size() == 0 || st.top() < lower || st.top() > upper) {
            return nullptr;
        }
        int val = st.top();
        st.pop();
        TreeNode::TreeNode *root = new TreeNode::TreeNode(val);
        root->right = construct(val, upper, st);
        root->left = construct(lower, val, st);
        return root;
    }

// Encodes a tree to a single string.
    string serialize(TreeNode::TreeNode *root) {
        string res;
        vector<int> arr;
        postOrder(root, arr);
        if (arr.size() == 0) {
            return res;
        }
        for (int i = 0; i < arr.size() - 1; i++) {
            res.append(to_string(arr[i]) + ",");
        }
        res.append(to_string(arr.back()));
        return res;
    }

// Decodes your encoded data to tree.
    TreeNode::TreeNode *deserialize(string data) {
        if (data.size() == 0) {
            return nullptr;
        }
        vector<string> arr = split(data, ',');
        stack<int> st;
        for (auto &str : arr) {
            st.emplace(stoi(str));
        }
        return construct(INT_MIN, INT_MAX, st);
    }
} // namespace SerializingAndDeserializingForBinaryTrees

void SerializingAndDeserializingForBinaryTrees_test() {
    string token = "213";
    TreeNode::TreeNode *root =
            SerializingAndDeserializingForBinaryTrees::deserialize(token);
    string ans = SerializingAndDeserializingForBinaryTrees::serialize(root);
    cout << ans << endl;
}

namespace deleteNode {
    TreeNode::TreeNode *dfs(TreeNode::TreeNode *root, int key) {
        if (root == nullptr)
            return nullptr;
        if (root->val == key) {
            if (root->left != nullptr) {
                root = root->left;
            } else if (root->right == nullptr) {
                root = root->right;
            } else {
                root = nullptr;
            }
            return root;
        } else {
            root = dfs(root->left, key);
            root = dfs(root->right, key);
        }
        return root;
    }

    TreeNode::TreeNode *deleteNode(TreeNode::TreeNode *root, int key) {
        if (root == nullptr) {
            return nullptr;
        }
        if (root->val > key) {
            root->left = deleteNode(root->left, key);
            return root;
        }
        if (root->val < key) {
            root->right = deleteNode(root->right, key);
            return root;
        }
        if (root->val == key) {
            if (!root->left && !root->right) {
                return nullptr;
            }
            if (!root->right) {
                return root->left;
            }
            if (!root->left) {
                return root->right;
            }
            TreeNode::TreeNode *successor = root->right;
            while (successor->left) {
                successor = successor->left;
            }
            root->right = deleteNode(root->right, successor->val);
            successor->right = root->right;
            successor->left = root->left;
            return successor;
        }
        return root;
    }
} // namespace deleteNode

void deleteNode_test() {
    vector<int> nums = {5, 3, 6, 2, 4, 0, 7};
    int key = 3;
    auto root = TreeNode::createTree(nums);
    auto ans = deleteNode::deleteNode(root, key);
    cout << TreeNode::print_tree(ans) << endl;
}

namespace frequencySort {
    bool cmp(pair<char, int> &p1, pair<char, int> &p2) {
        return p1.second < p2.second;
    }

    string frequencySort(string s) {
        string ans;
        unordered_map<char, int> umap;
        for (auto &c : s) {
            umap[c]++;
        }
        vector<pair<char, int>> map;
        for (auto p : umap) {
            map.push_back(p);
        }

        sort(map.begin(), map.end(), [](pair<char, int> &p1, pair<char, int> &p2) {
            return p1.second > p2.second;
        });
        for (auto tmp : map) {
            int index = 0;
            while (index < tmp.second) {
                index++;
                ans += tmp.first;
            }
        }
        return ans;
    }
} // namespace frequencySort
void frequencySort_test() {
    string s = "2a554442f544asfasssffffasss";
    cout << s << " 根据字符出现频率排序后 " << frequencySort::frequencySort(s)
         << endl;
}

namespace minMoves {
    int minMoves(vector<int> &nums) {
        int min_num = *min_element(nums.begin(), nums.end());
        int ans = 0;
        for (auto num : nums) {
            ans += num - min_num;
        }
        return ans;
    }
} // namespace minMoves

void minMoves_test() {
    vector<int> nums = {1, 2, 3};
    cout << "最小操作次数：" << minMoves::minMoves(nums) << endl;
}

namespace fourSumCount {
    int fourSumCount(vector<int> &nums1, vector<int> &nums2, vector<int> &nums3,
                     vector<int> &nums4) {
        unordered_map<int, int> count12;
        for (int u : nums1) {
            for (int v : nums2) {
                ++count12[u + v];
            }
        }
        int ans = 0;
        for (int u : nums3) {
            for (int v : nums4) {
                if (count12.count(-u - v)) {
                    ans += count12[-u - v];
                }
            }
        }
        return ans;
    }
} // namespace fourSumCount

void fourSumCount_test() {
    vector<int> nums1 = {1, 2}, nums2 = {-2, -1}, nums3 = {-1, 2}, nums4 = {0, 2};
    cout << "四数相加为0的元组数："
         << fourSumCount::fourSumCount(nums1, nums2, nums3, nums4) << endl;
}

namespace findContentChildren {
    int findContentChildren(vector<int> &g, vector<int> &s) {
        sort(g.begin(), g.end());
        sort(s.begin(), s.end());
        int m = g.size(), n = s.size();
        int count = 0;
        for (int i = 0, j = 0; i < m && j < n; i++, j++) {
            while (j < n && g[i] > s[j]) {
                j++;
            }
            if (j < n) {
                count++;
            }
        }
        return count;
    }
} // namespace findContentChildren

void findContentChildren_test() {
    vector<int> g = {1, 2, 3};
    vector<int> s = {1, 1};
    cout << "有 " << findContentChildren::findContentChildren(g, s)
         << " 个小孩儿被满足" << endl;
}

namespace find132pattern {
    bool find132pattern(vector<int> &nums) {
        int n = nums.size();
        if (n < 3) {
            return false;
        }

        // 左侧最小值
        int left_min = nums[0];
        // 右侧所有元素
        multiset<int> right_all;

        for (int k = 2; k < n; ++k) {
            right_all.insert(nums[k]);
        }

        for (int j = 1; j < n - 1; ++j) {
            if (left_min < nums[j]) {
                auto it = right_all.upper_bound(left_min);
                if (it != right_all.end() && *it < nums[j]) {
                    return true;
                }
            }
            left_min = min(left_min, nums[j]);
            right_all.erase(right_all.find(nums[j + 1]));
        }

        return false;
    }
} // namespace find132pattern

void find132pattern_test() {
    vector<int> nums = {3, 5, 0, 3, 4};
    string s = "是";
    if (find132pattern::find132pattern(nums)) {
        s = "存在";
    } else
        s = "不存在";
    cout << "nums " << s << " 132模式子序列" << endl;
}

namespace hammingDistance {
    int hammingDistance(int x, int y) {
        int xory = x ^y;
        int ans = 0;
        while (xory != 0) {
            ans += xory & 1;
            xory >>= 1;
        }
        return ans;
    }
} // namespace hammingDistance

void hammingDistance_test() {
    int x = 1;
    int y = 4;
    cout << x << "和" << y << "的汉明距离："
         << hammingDistance::hammingDistance(x, y) << endl;
}

namespace minMoves2 {
    int partition(std::vector<int> &nums, int left, int right) {
        int pivotIndex = left + rand() % (right - left + 1);
        int pivotValue = nums[pivotIndex];
        std::swap(nums[pivotIndex], nums[right]);
        int storeIndex = left;

        for (int i = left; i < right; ++i) {
            if (nums[i] < pivotValue) {
                std::swap(nums[i], nums[storeIndex]);
                ++storeIndex;
            }
        }
        std::swap(nums[storeIndex], nums[right]);
        return storeIndex;
    }

    int quickSelect(std::vector<int> &nums, int left, int right, int k) {
        if (left == right) {
            return nums[left];
        }

        int pivotIndex = partition(nums, left, right);

        if (k == pivotIndex) {
            return nums[k];
        } else if (k < pivotIndex) {
            return quickSelect(nums, left, pivotIndex - 1, k);
        } else {
            return quickSelect(nums, pivotIndex + 1, right, k);
        }
    }

    int findKthSmallest(std::vector<int> &nums, int k) {
        return quickSelect(nums, 0, nums.size() - 1, k - 1);
    }

    int minMoves2(vector<int> &nums) {
        sort(nums.begin(), nums.end());

        int n = nums.size();
        int x = findKthSmallest(nums, n / 2);
        int ret = 0;
        for (int i = 0; i < n; i++) {
            ret += abs(nums[i] - x);
        }
        return ret;
    }
} // namespace minMoves2

void minMoves2_test() {
    vector<int> nums = {1, 10, 2, 9};
    cout << "nums最小移动距离2：" << minMoves2::minMoves2(nums) << endl;
}

namespace islandPerimeter {
    constexpr static int dx[4] = {0, 1, 0, -1};
    constexpr static int dy[4] = {1, 0, -1, 0};

    int dfs(int x, int y, vector<vector<int>> &grid, int n, int m) {
        if (x < 0 || x >= n || y < 0 || y >= m || grid[x][y] == 0) {
            return 1;
        }
        if (grid[x][y] == 2) {
            return 0;
        }
        grid[x][y] = 2;
        int res = 0;
        for (int i = 0; i < 4; ++i) {
            int tx = x + dx[i];
            int ty = y + dy[i];
            res += dfs(tx, ty, grid, n, m);
        }
        return res;
    }

    int islandPerimeter(vector<vector<int>> &grid) {
        int n = grid.size(), m = grid[0].size();
        int ans = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (grid[i][j] == 1) {
                    ans += dfs(i, j, grid, n, m);
                }
            }
        }
        return ans;
    }
} // namespace islandPerimeter

void islandPerimeter_test() {
    vector<vector<int>> grid = {
            {0, 1, 0, 0},
            {1, 1, 1, 0},
            {0, 1, 0, 0},
            {1, 1, 0, 0}};
    cout << "岛屿grid的边界长度：" << islandPerimeter::islandPerimeter(grid)
         << endl;
}

namespace canIWin {
    bool canIWinHelper(int maxChoosableInteger, int desiredTotal, int chosen,
                       std::unordered_map<int, bool> &memo) {
        if (memo.count(chosen)) {
            return memo[chosen];
        }

        for (int i = 1; i <= maxChoosableInteger; ++i) {
            int mask = 1 << (i - 1);
            if ((chosen & mask) == 0) {
                if (i >= desiredTotal ||
                    !canIWinHelper(maxChoosableInteger, desiredTotal - i, chosen | mask,
                                   memo)) {
                    memo[chosen] = true;
                    return true;
                }
            }
        }

        memo[chosen] = false;
        return false;
    }

    bool canIWin(int maxChoosableInteger, int desiredTotal) {
        if (maxChoosableInteger >= desiredTotal) {
            return true;
        }

        int sum = (1 + maxChoosableInteger) * maxChoosableInteger / 2;
        if (sum < desiredTotal) {
            return false;
        }

        std::unordered_map<int, bool> memo;
        return canIWinHelper(maxChoosableInteger, desiredTotal, 0, memo);
    }
} // namespace canIWin

void canIWin_test() {
    int maxChoosableInteger = 10;
    int desiredTotal = 21;
    string ans = "是";
    if (canIWin::canIWin(maxChoosableInteger, desiredTotal))
        ans = "必胜";
    else
        ans = "必败";
    cout << "可选取maxChoosableInteger：" << maxChoosableInteger << "，目标数值："
         << desiredTotal << ",先手方 " << ans << endl;
}

namespace findSubstringInWraproundString {
    int findSubstringInWraproundString(string s) {
        vector<int> dp(26);
        int k = 0;
        for (int i = 0; i < s.length(); ++i) {
            if (i && (s[i] - s[i - 1] + 26) % 26 == 1) { // 字符之差为 1 或 -25
                ++k;
            } else {
                k = 1;
            }
            dp[s[i] - 'a'] = max(dp[s[i] - 'a'], k);
        }
        return accumulate(dp.begin(), dp.end(), 0);
    }
} // namespace findSubstringInWraproundString

void findSubstringInWraproundString_test() {
    string s = "zab";
    cout << "字符串" << s << " 有 "
         << findSubstringInWraproundString::findSubstringInWraproundString(s)
         << " 个不同子串" << endl;
}

namespace findAllConcatenatedWordsInADict {
    struct Trie {
        bool isEnd;
        vector<Trie *> children;

        Trie() {
            this->children = vector<Trie *>(26, nullptr);
            this->isEnd = false;
        }
    };

    Trie *trie = new Trie();

    bool dfs(const string &word, int start, vector<int> &visited) {
        if (word.size() == start) {
            return true;
        }
        if (visited[start]) {
            return false;
        }
        visited[start] = true;
        Trie *node = trie;
        for (int i = start; i < word.size(); i++) {
            char ch = word[i];
            int index = ch - 'a';
            node = node->children[index];
            if (node == nullptr) {
                return false;
            }
            if (node->isEnd) {
                if (dfs(word, i + 1, visited)) {
                    return true;
                }
            }
        }
        return false;
    }

    void insert(const string &word) {
        Trie *node = trie;
        for (int i = 0; i < word.size(); i++) {
            char ch = word[i];
            int index = ch - 'a';
            if (node->children[index] == nullptr) {
                node->children[index] = new Trie();
            }
            node = node->children[index];
        }
        node->isEnd = true;
    }

    vector<string> findAllConcatenatedWordsInADict(vector<string> &words) {
        vector<string> ans;
        sort(words.begin(), words.end(),
             [&](const string &a, const string &b) { return a.size() < b.size(); });
        for (int i = 0; i < words.size(); i++) {
            string word = words[i];
            if (word.size() == 0) {
                continue;
            }
            vector<int> visited(word.size(), 0);
            if (dfs(word, 0, visited)) {
                ans.emplace_back(word);
            } else {
                insert(word);
            }
        }
        return ans;
    }
} // namespace findAllConcatenatedWordsInADict

void findAllConcatenatedWordsInADict_test() {
    vector<string> words = {"cat", "cats", "catsdogcats",
                            "dog", "dogcatsdog", "hippopotamuses",
                            "rat", "ratcatdogcat"};
    auto ans =
            findAllConcatenatedWordsInADict::findAllConcatenatedWordsInADict(words);

    for (auto word : ans) {
        cout << word << " ";
    }
    cout << endl;
}

namespace makesquare {
    bool dfs(int index, vector<int> &matchsticks, vector<int> &edges, int len) {
        if (index == matchsticks.size()) {
            return true;
        }
        for (int i = 0; i < edges.size(); i++) {
            edges[i] += matchsticks[index];
            if (edges[i] <= len && dfs(index + 1, matchsticks, edges, len)) {
                return true;
            }
            edges[i] -= matchsticks[index];
        }
        return false;
    }

    bool makesquare(vector<int> &matchsticks) {
        int totalLen = accumulate(matchsticks.begin(), matchsticks.end(), 0);
        if (totalLen % 4 != 0) {
            return false;
        }
        sort(matchsticks.begin(), matchsticks.end(), greater<int>()); // 减少搜索量

        vector<int> edges(4);
        return dfs(0, matchsticks, edges, totalLen / 4);
    }
} // namespace makesquare

void makesquare_test() {
    vector<int> matchsticks = {3, 3, 3, 3, 4};
    cout << makesquare::makesquare(matchsticks) << endl;
}

namespace findMaxForm {
    vector<int> getZerosOnes(string &str) {
        vector<int> zeroOnes(2);
        int length = str.length();
        for (int i = 0; i < length; ++i) {
            zeroOnes[str[i] - '0']++;
        }
        return zeroOnes;
    }

    int findMaxForm(vector<string> &strs, int m, int n) {
        int length = strs.size();
        vector<vector<vector<int >>> dp(
                length + 1, vector<vector<int >>(m + 1, vector<int>(n + 1)));
        for (int i = 1; i <= length; i++) {
            vector<int> &&zerosOnes = getZerosOnes(strs[i - 1]);
            int zeros = zerosOnes[0], ones = zerosOnes[1];
            for (int j = 0; j <= m; ++j) {
                for (int k = 0; k <= n; ++k) {
                    dp[i][j][k] = dp[i - 1][j][k];
                    if (j >= zeros && k >= ones) {
                        dp[i][j][k] = max(dp[i][j][k], dp[i - 1][j - zeros][k - ones] + 1);
                    }
                }
            }
        }
        return dp[length][m][n];
    }
} // namespace findMaxForm

void findMaxForm_test() {
    vector<string> strs = {"10", "0001", "111001", "1", "0"};
    int m = 5, n = 3;
    cout << findMaxForm::findMaxForm(strs, m, n) << endl;
}

namespace findRadius {
    int findRadius(vector<int> &houses, vector<int> &heaters) {
        sort(houses.begin(), houses.end());
        sort(heaters.begin(), heaters.end());
        int ans = 0;
        for (int i = 0, j = 0; i < houses.size(); i++) {
            int curDistance = abs(houses[i] - heaters[j]);
            while (j < heaters.size() - 1 &&
                   abs(houses[i] - heaters[j]) >= abs(houses[i] - heaters[j + 1])) {
                j++;
                curDistance = min(curDistance, abs(houses[i] - heaters[j]));
            }
            ans = max(ans, curDistance);
        }
        return ans;
    }
} // namespace findRadius

void findRadius_test() {
    vector<int> houses = {1, 2, 3};
    vector<int> heaters = {2};
    // cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
    houses = {1, 2, 3, 4};
    heaters = {1, 4};
    // cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
    houses = {1, 5};
    heaters = {2};
    // cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
    houses = {1, 5};
    heaters = {10};
    cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
}

#include <ctime>

namespace randomlyGeneratePointsWithinACircle {

    class Solution {
        mt19937 gen{random_device{}()};
        uniform_real_distribution<double> dis;
        double xc, yc, r;

    public:
        Solution(double radius, double x_center, double y_center)
                : dis(-radius, radius), xc(x_center), yc(y_center), r(radius) {}

        vector<double> randPoint() {
            while (true) {
                double x = dis(gen), y = dis(gen);
                if (x * x + y * y <= r * r) {
                    return {xc + x, yc + y};
                }
            }
        }
    };
} // namespace randomlyGeneratePointsWithinACircle

void randomlyGeneratePointsWithinACircle_test() {
    vector<double> param = {1.0, 0.0, 0.0};
    randomlyGeneratePointsWithinACircle::Solution sol(param[0], param[1],
                                                      param[2]);
    cout << sol.randPoint()[0] << "," << sol.randPoint()[1] << endl;
    cout << sol.randPoint()[0] << "," << sol.randPoint()[1] << endl;
    cout << sol.randPoint()[0] << "," << sol.randPoint()[1] << endl;
}

namespace largestPalindrome {
    int largestPalindrome(int n) {
        if (n == 1) {
            return 9;
        }
        int upper = pow(10, n) - 1;
        for (int left = upper;; --left) { // 枚举回文数的左半部分
            long p = left;
            for (int x = left; x > 0; x /= 10) {
                p = p * 10 + x % 10; // 翻转左半部分到其自身末尾，构造回文数 p
            }
            for (long x = upper; x * x >= p; --x) {
                if (p % x == 0) { // x 是 p 的因子
                    return p % 1337;
                }
            }
        }
    }
} // namespace largestPalindrome

void largestPalindrome_test() {
    int n = 2;
    cout << largestPalindrome::largestPalindrome(2) << endl;
}

namespace medianSlidingWindow {
    class DualHeap {
    private:
        // 大根堆，维护较小的一半元素
        priority_queue<int> small;
        // 小根堆，维护较大的一半元素
        priority_queue<int, vector<int>, greater<int>> large;
        // 哈希表，记录「延迟删除」的元素，key 为元素，value 为需要删除的次数
        unordered_map<int, int> delayed;

        int k;
        // small 和 large 当前包含的元素个数，需要扣除被「延迟删除」的元素
        int smallSize, largeSize;

    public:
        DualHeap(int _k) : k(_k), smallSize(0), largeSize(0) {}

    private:
        // 不断地弹出 heap 的堆顶元素，并且更新哈希表
        template<typename T>
        void prune(T &heap) {
            while (!heap.empty()) {
                int num = heap.top();
                if (delayed.count(num)) {
                    --delayed[num];
                    if (!delayed[num]) {
                        delayed.erase(num);
                    }
                    heap.pop();
                } else {
                    break;
                }
            }
        }

        // 调整 small 和 large 中的元素个数，使得二者的元素个数满足要求
        void makeBalance() {
            if (smallSize > largeSize + 1) {
                // small 比 large 元素多 2 个
                large.push(small.top());
                small.pop();
                --smallSize;
                ++largeSize;
                // small 堆顶元素被移除，需要进行 prune
                prune(small);
            } else if (smallSize < largeSize) {
                // large 比 small 元素多 1 个
                small.push(large.top());
                large.pop();
                ++smallSize;
                --largeSize;
                // large 堆顶元素被移除，需要进行 prune
                prune(large);
            }
        }

    public:
        void insert(int num) {
            if (small.empty() || num <= small.top()) {
                small.push(num);
                ++smallSize;
            } else {
                large.push(num);
                ++largeSize;
            }
            makeBalance();
        }

        void erase(int num) {
            ++delayed[num];
            if (num <= small.top()) {
                --smallSize;
                if (num == small.top()) {
                    prune(small);
                }
            } else {
                --largeSize;
                if (num == large.top()) {
                    prune(large);
                }
            }
            makeBalance();
        }

        double getMedian() {
            return k & 1 ? small.top() : ((double) small.top() + large.top()) / 2;
        }
    };

    vector<double> medianSlidingWindow(vector<int> &nums, int k) {
        DualHeap dh(k);
        for (int i = 0; i < k; ++i) {
            dh.insert(nums[i]);
        }
        vector<double> ans = {dh.getMedian()};
        for (int i = k; i < nums.size(); ++i) {
            dh.insert(nums[i]);
            dh.erase(nums[i - k]);
            ans.push_back(dh.getMedian());
        }
        return ans;
    }
} // namespace medianSlidingWindow

void medianSlidingWindow_test() {
    vector<int> nums = {1, 3, -1, -3, 5, 3, 6, 7};
    int k = 3;
    auto ans = medianSlidingWindow::medianSlidingWindow(nums, k);
    print_vector(ans);
    nums = {1, 2, 3, 4, 2, 3, 1, 4, 2};
    ans = medianSlidingWindow::medianSlidingWindow(nums, k);
    print_vector(ans);
}

namespace findMaxConsecutiveOnes {
    int findMaxConsecutiveOnes(vector<int> &nums) {
        int ans = 0;
        int tmp = 0;
        for (auto num : nums) {
            if (num == 1) {
                tmp++;
                if (tmp > ans)
                    ans = tmp;
            } else {
                tmp = 0;
            }
        }
        return ans;
    }
} // namespace findMaxConsecutiveOnes

void findMaxConsecutiveOnes_test() {
    vector<int> nums = {1, 1, 0, 1, 1, 1};
    cout << findMaxConsecutiveOnes::findMaxConsecutiveOnes(nums) << endl;
    nums = {1, 0, 1, 1, 0, 1};
    cout << findMaxConsecutiveOnes::findMaxConsecutiveOnes(nums) << endl;
}

namespace predictTheWinner {
    int total(vector<int> &nums, int start, int end, int turn) {
        if (start == end) {
            return nums[start] * turn;
        }
        int scoreStart = nums[start] * turn + total(nums, start + 1, end, -turn);
        int scoreEnd = nums[end] * turn + total(nums, start, end - 1, -turn);
        return max(scoreStart * turn, scoreEnd * turn) * turn;
    }

    bool predictTheWinner(vector<int> &nums) {
        return total(nums, 0, nums.size() - 1, 1) >= 0;
    }
} // namespace predictTheWinner

void predictTheWinner_test() {
    vector<int> nums = {1, 5, 2};
    cout << predictTheWinner::predictTheWinner(nums) << endl;
    nums = {1, 5, 233, 7};
    cout << predictTheWinner::predictTheWinner(nums) << endl;
}

namespace findMinStep {
    struct State {
        string board;
        string hand;
        int step;

        State(const string &board, const string &hand, int step) {
            this->board = board;
            this->hand = hand;
            this->step = step;
        }
    };

    string clean(const string &s) {
        string res;
        vector<pair<char, int>> st;

        for (auto c : s) {
            while (!st.empty() && c != st.back().first && st.back().second >= 3) {
                st.pop_back();
            }
            if (st.empty() || c != st.back().first) {
                st.push_back({c, 1});
            } else {
                st.back().second++;
            }
        }
        if (!st.empty() && st.back().second >= 3) {
            st.pop_back();
        }
        for (int i = 0; i < st.size(); ++i) {
            for (int j = 0; j < st[i].second; ++j) {
                res.push_back(st[i].first);
            }
        }
        return res;
    }

    int findMinStep(string board, string hand) {
        unordered_set<string> visited;
        sort(hand.begin(), hand.end());

        visited.insert(board + " " + hand);
        queue<State> qu;
        qu.push(State(board, hand, 0));
        while (!qu.empty()) {
            State curr = qu.front();
            qu.pop();

            for (int j = 0; j < curr.hand.size(); ++j) {
                // 第 1 个剪枝条件: 当前选择的球的颜色和前一个球的颜色相同
                if (j > 0 && curr.hand[j] == curr.hand[j - 1]) {
                    continue;
                }
                for (int i = 0; i <= curr.board.size(); ++i) {
                    // 第 2 个剪枝条件: 只在连续相同颜色的球的开头位置插入新球
                    if (i > 0 && curr.board[i - 1] == curr.hand[j]) {
                        continue;
                    }

                    // 第 3 个剪枝条件: 只在以下两种情况放置新球
                    bool choose = false;
                    //   第 1 种情况 : 当前球颜色与后面的球的颜色相同
                    if (i < curr.board.size() && curr.board[i] == curr.hand[j]) {
                        choose = true;
                    }
                    //   第 2 种情况 : 当前后颜色相同且与当前颜色不同时候放置球
                    if (i > 0 && i < curr.board.size() &&
                        curr.board[i - 1] == curr.board[i] &&
                        curr.board[i] != curr.hand[j]) {
                        choose = true;
                    }
                    if (choose) {
                        string new_board = clean(curr.board.substr(0, i) + curr.hand[j] +
                                                 curr.board.substr(i));
                        string new_hand = curr.hand.substr(0, j) + curr.hand.substr(j + 1);
                        if (new_board.size() == 0) {
                            return curr.step + 1;
                        }
                        if (!visited.count(new_board + " " + new_hand)) {
                            qu.push(State(new_board, new_hand, curr.step + 1));
                            visited.insert(new_board + " " + new_hand);
                        }
                    }
                }
            }
        }

        return -1;
    }
} // namespace findMinStep

void findMinStep_test() {
    string board = "WRRBBW";
    string hand = "RB";
    cout << board << "," << hand << "," << findMinStep::findMinStep(board, hand)
         << endl;
    board = "WWRRBBWW";
    hand = "WRBRW";
    cout << board << "," << hand << "," << findMinStep::findMinStep(board, hand)
         << endl;
    board = "G";
    hand = "GGGGG";
    cout << board << "," << hand << "," << findMinStep::findMinStep(board, hand)
         << endl;
}

namespace nextGreaterElement {
    vector<int> nextGreaterElement(vector<int> &nums1, vector<int> &nums2) {
        unordered_map<int, int> map;
        stack<int> st;
        for (int i = nums2.size() - 1; i >= 0; --i) {
            auto num = nums2[i];
            while (!st.empty() && num >= st.top()) {
                st.pop();
            }
            map[num] = st.empty() ? -1 : st.top();
            st.push(num);
        }
        vector<int> res(nums1.size());
        for (int i = 0; i < nums1.size(); ++i) {
            res[i] = map[nums1[i]];
        }
        return res;
    }
} // namespace nextGreaterElement

void nextGreaterElement_test() {
    vector<int> nums1 = {4, 1, 2};
    vector<int> nums2 = {1, 3, 4, 2};
    auto res = nextGreaterElement::nextGreaterElement(nums1, nums2);
    print_vector(res);
    nums1 = {1, 3, 5, 2, 4};;
    nums2 = {6, 5, 4, 3, 2, 1, 7};
    res = nextGreaterElement::nextGreaterElement(nums1, nums2);
    print_vector(res);
}

namespace findDiagonalOrder {
    vector<int> findDiagonalOrder(vector<vector<int>> &mat) {
        int m = mat.size();
        int n = mat[0].size();
        vector<int> res;
        for (int i = 0; i < m + n - 1; i++) {
            if (i % 2) {
                int x = i < n ? 0 : i - n + 1;
                int y = i < n ? i : n - 1;
                while (x < m && y >= 0) {
                    res.emplace_back(mat[x][y]);
                    x++;
                    y--;
                }
            } else {
                int x = i < m ? i : m - 1;
                int y = i < m ? 0 : i - m + 1;
                while (x >= 0 && y < n) {
                    res.emplace_back(mat[x][y]);
                    x--;
                    y++;
                }
            }
        }
        return res;
    }
} // namespace findDiagonalOrder

void findDiagonalOrder_test() {
    vector<vector<int>> mat = {{1, 2, 3},
                               {4, 5, 6},
                               {7, 8, 9}};
    auto ans = findDiagonalOrder::findDiagonalOrder(mat);
    print_vector(ans);
}

namespace findWords {
    vector<string> findWords(vector<string> &words) {
        vector<string> ans;
        unordered_map<char, int> map;
        string line1 = "qwertyuiopQWERTYUIOP";
        string line2 = "asdfghjklASDFGHJKL";
        string line3 = "zxcvbnmZXCVBNM";
        for (auto c : line1) {
            map[c] = 1;
        }
        for (auto c : line2) {
            map[c] = 2;
        }
        for (auto c : line3) {
            map[c] = 3;
        }
        for (auto s : words) {
            string tmp;
            int flag = map[s[0]];
            for (auto c : s) {
                if (map[c] != flag)
                    break;
                else
                    tmp += c;
            }
            if (tmp.size() == s.size())
                ans.push_back(tmp);
        }
        return ans;
    }
} // namespace findWords

void findWords_test() {
    vector<string> words, ans;
    words = {"Hello", "Alaska", "Dad", "Peace"};
    ans = findWords::findWords(words);
    print_vector(ans);
    words = {"omk"};
    ans = findWords::findWords(words);
    print_vector(ans);
    words = {"adsdf", "sfd"};
    ans = findWords::findWords(words);
    print_vector(ans);
}

namespace findMode {
    vector<int> answer;
    int base, count, maxCount;

    void update(int x) {
        if (x == base) {
            ++count;
        } else {
            count = 1;
            base = x;
        }
        if (count == maxCount) {
            answer.push_back(base);
        }
        if (count > maxCount) {
            maxCount = count;
            answer = vector<int>{base};
        }
    }

    void dfs(TreeNode::TreeNode *root) {
        if (root == nullptr)
            return;
        dfs(root->left);
        update(root->val);
        dfs(root->right);
    }

    vector<int> findMode(TreeNode::TreeNode *root) {
        answer.clear();
        base = 0;
        count = 0;
        maxCount = 0;
        dfs(root);
        return answer;
    }
} // namespace findMode

void findMode_test() {
    vector<int> tree = {1, 0, 2, 2};
    TreeNode::TreeNode *root;
    root = create_treenode(tree);
    vector<int> ans;
    ans = findMode::findMode(root);
    print_vector(ans);
    tree = {0};
    root = create_treenode(tree);
    ans = findMode::findMode(root);
    print_vector(ans);
}

namespace findMaximizedCapital {
    typedef pair<int, int> pii;

    int findMaximizedCapital(int k, int w, vector<int> &profits,
                             vector<int> &capital) {
        int n = profits.size();
        int curr = 0;
        priority_queue<int, vector<int>, less<int>> pq;
        vector<pii> arr;
        for (int i = 0; i < n; ++i) {
            arr.push_back({capital[i], profits[i]});
        }
        sort(arr.begin(), arr.end());
        for (int i = 0; i < k; ++i) {
            while (curr < n && arr[curr].first <= w) {
                pq.push(arr[curr].second);
                curr++;
            }
            if (!pq.empty()) {
                w += pq.top();
                pq.pop();
            } else {
                break;
            }
        }
        return w;
    }
} // namespace findMaximizedCapital

void findMaximizedCapital_test() {
    int k, w;
    k = 2;
    w = 0;
    vector<int> profits, captial;
    profits = {1, 2, 3};
    captial = {0, 1, 1};
    cout << findMaximizedCapital::findMaximizedCapital(k, w, profits, captial)
         << endl;
    k = 3;
    w = 0;
    profits = {1, 2, 3};
    captial = {0, 1, 2};
    cout << findMaximizedCapital::findMaximizedCapital(k, w, profits, captial)
         << endl;
}

namespace nextGreaterElements {
    vector<int> nextGreaterElements(vector<int> &nums) {
        int n = nums.size();
        stack<int> stk;
        vector<int> ret(n, -1);
        for (int i = 0; i < 2 * n - 1; ++i) {
            while (!stk.empty() && nums[stk.top()] < nums[i % n]) {
                ret[stk.top()] = nums[i % n];
                stk.pop();
            }
            stk.push(i % n);
        }
        return ret;
    }
} // namespace nextGreaterElements

void nextGreaterElements_test() {
    vector<int> ans, nums;
    nums = {1, 2, 1};
    ans = nextGreaterElements::nextGreaterElements(nums);
    print_vector(ans);
    nums = {1, 2, 3, 4, 3};
    ans = nextGreaterElements::nextGreaterElements(nums);
    print_vector(ans);
    nums = {100, 1, 11, 1, 120, 111, 123, 1, -1, -100};
    ans = nextGreaterElements::nextGreaterElements(nums);
    print_vector(ans);
}

namespace convertToBase7 {
    string convertToBase7(int num) {
        if (num == 0) {
            return "0";
        }
        bool negative = num < 0;
        num = abs(num);
        string digits;
        while (num > 0) {
            digits.push_back(num % 7 + '0');
            num /= 7;
        }
        if (negative) {
            digits.push_back('-');
        }
        reverse(digits.begin(), digits.end());
        return digits;
    }
} // namespace convertToBase7

void convertToBase7_test() {
    int num;
    num = 100;
    cout << num << " 的7进制数为 " << convertToBase7::convertToBase7(num) << endl;
    num = -7;
    cout << num << " 的7进制数为 " << convertToBase7::convertToBase7(num) << endl;
}

namespace findRelativeRanks {
    typedef pair<int, int> pii;

    vector<string> findRelativeRanks(vector<int> &score) {
        vector<pii> map;
        vector<string> ans(score.size());
        for (int i = 0; i < score.size(); ++i) {
            map.push_back({score[i], i});
        }
        sort(map.begin(), map.end(), greater<pii>());

        for (int i = 0; i < map.size(); ++i) {
            auto tmp = map[i];
            if (i == 0) {
                ans[tmp.second] = "Gold Medal";
            } else if (i == 1) {
                ans[tmp.second] = "Silver Medal";
            } else if (i == 2) {
                ans[tmp.second] = "Bronze Medal";
            } else {
                ans[tmp.second] = to_string(i + 1);
            }
        }
        return ans;
    }
} // namespace findRelativeRanks

void findRelativeRanks_test() {
    vector<int> score;
    vector<string> ans;
    score = {5, 4, 3, 2, 1};
    ans = findRelativeRanks::findRelativeRanks(score);
    print_vector(ans);
    score = {10, 3, 8, 9, 4};
    ans = findRelativeRanks::findRelativeRanks(score);
    print_vector(ans);
}

namespace findFrequentTreeSum {
    unordered_map<int, int> cnt;
    int maxCnt = 0;

    int dfs(TreeNode::TreeNode *node) {
        if (node == nullptr) {
            return 0;
        }
        int sum = node->val + dfs(node->left) + dfs(node->right);
        maxCnt = max(maxCnt, ++cnt[sum]);
        return sum;
    }

    vector<int> findFrequentTreeSum(TreeNode::TreeNode *root) {
        dfs(root);
        vector<int> ans;
        for (auto &[s, c] : cnt) {
            if (c == maxCnt) {
                ans.emplace_back(s);
            }
        }
        return ans;
    }
} // namespace findFrequentTreeSum

void findFrequentTreeSum_test() {
    TreeNode::TreeNode *root = create_treenode({5, 2, -3});
    vector<int> ans;
    ans = findFrequentTreeSum::findFrequentTreeSum(root);
    print_vector(ans);
    root = create_treenode({5, 2, -5});
    ans = findFrequentTreeSum::findFrequentTreeSum(root);
    print_vector(ans);
}

namespace findBottomLeftValue {
    int findBottomLeftValue(TreeNode::TreeNode *root) {
        if (root == nullptr) {
            return 0;
        }
        // bfs
        queue<TreeNode::TreeNode *> stk;
        stk.push(root);
        int size = stk.size();
        unordered_map<int, int> map;
        int laynum = 0;
        int ans = root->val;
        while (stk.size() != 0) {
            size = stk.size();
            for (int i = 0; i < size; ++i) {
                TreeNode::TreeNode *p = stk.front();
                stk.pop();
                map[laynum] = p->val;
                if (i == 0)
                    ans = p->val;
                if (p->left)
                    stk.push(p->left);
                if (p->right)
                    stk.push(p->right);
            }
            laynum++;
        }
        return ans;
    }
} // namespace findBottomLeftValue

void findBottomLeftValue_test() {
    TreeNode::TreeNode *root;
    root = create_treenode({2, 1, 3});
    cout << findBottomLeftValue::findBottomLeftValue(root) << endl;
    root = create_treenode({1, 2, 3, 4, 0, 5, 6, 0, 0, 7});
    cout << findBottomLeftValue::findBottomLeftValue(root) << endl;
}

namespace longestPalindromeSubseq {
    int longestPalindromeSubseq(string s) {
        int n = s.size();
        vector<vector<int>> dp(n, vector<int>(n));
        for (int i = n - 1; i >= 0; --i) {
            char c1 = s[i];
            dp[i][i] = 1;
            for (int j = i + 1; j < n; ++j) {
                char c2 = s[j];
                if (c1 == c2) {
                    dp[i][j] = dp[i + 1][j - 1] + 2;
                } else {
                    dp[i][j] = max(dp[i + 1][j], dp[i][j - 1]);
                }
            }
        }
        return dp[0][n - 1];
    }
} // namespace longestPalindromeSubseq

void longestPalindromeSubseq_test() {
    string s = "bbbab";
    cout << s << " 的最大回文子串序列个数："
         << longestPalindromeSubseq::longestPalindromeSubseq(s) << endl;
    s = "cbbd";
    cout << s << " 的最大回文子串序列个数："
         << longestPalindromeSubseq::longestPalindromeSubseq(s) << endl;
}

namespace findMinMoves {
    int findMinMoves(vector<int> &machines) {
        int tot = accumulate(machines.begin(), machines.end(), 0);
        int n = machines.size();
        if (tot % n) {
            return -1;
        }
        int avg = tot / n;
        int ans = 0, sum = 0;
        for (auto num : machines) {
            num -= avg;
            sum += num;
            ans = max(ans, max(abs(sum), num));
        }
        return ans;
    }
} // namespace findMinMoves

void findMinMoves_test() {
    vector<int> machines;
    machines = {1, 0, 5};
    cout << findMinMoves::findMinMoves(machines) << endl;
    machines = {0, 3, 0};
    cout << findMinMoves::findMinMoves(machines) << endl;
    machines = {0, 2, 0};
    cout << findMinMoves::findMinMoves(machines) << endl;
}

namespace change {
    int change(int amount, vector<int> &coins) {
        vector<int> dp(amount + 1);
        int n = coins.size();
        dp[0] = 1;
        for (auto &coin : coins) {
            for (int i = coin; i <= amount; ++i) {
                dp[i] += dp[i - coin];
            }
        }
        return dp[amount];
    }
} // namespace change

void change_test() {
    int amount;
    vector<int> coins;
    amount = 5;
    coins = {1, 2, 5};
    cout << "用coins ";
    print_vector(coins);
    cout << " 可以有：" << change::change(amount, coins) << " 种方式" << endl;
    amount = 3;
    coins = {2};
    cout << "用coins ";
    print_vector(coins);
    cout << " 可以有：" << change::change(amount, coins) << " 种方式" << endl;
    amount = 10;
    coins = {10};
    cout << "用coins ";
    print_vector(coins);
    cout << " 可以有：" << change::change(amount, coins) << " 种方式" << endl;
}

namespace findLUSlength {
    int findLUSlength(string a, string b) {
        return a != b ? max(a.length(), b.length()) : -1;
    }
} // namespace findLUSlength

void findLUSlength_test() {
    string a, b;
    a = "aba", b = "cdc";
    cout << "序列a " << a << " 和序列b " << b << " 两个字符串的最长特殊序列个数："
         << findLUSlength::findLUSlength(a, b) << endl;
    a = "aaa", b = "bbb";
    cout << "序列a " << a << " 和序列b " << b << " 两个字符串的最长特殊序列个数："
         << findLUSlength::findLUSlength(a, b) << endl;
    a = "aaa", b = "aaa";
    cout << "序列a " << a << " 和序列b " << b << " 两个字符串的最长特殊序列个数："
         << findLUSlength::findLUSlength(a, b) << endl;
}

namespace findLUSlength2 {
    int findLUSlength(vector<string> &strs) {
        auto is_subseq = [](const string &s, const string &t) -> bool {
            int pt_s = 0, pt_t = 0;
            while (pt_s < s.size() && pt_t < t.size()) {
                if (s[pt_s] == t[pt_t]) {
                    ++pt_s;
                }
                ++pt_t;
            }
            return pt_s == s.size();
        };

        int n = strs.size();
        int ans = -1;
        for (int i = 0; i < n; ++i) {
            bool check = true;
            for (int j = 0; j < n; ++j) {
                if (i != j && is_subseq(strs[i], strs[j])) {
                    check = false;
                    break;
                }
            }
            if (check) {
                ans = max(ans, static_cast<int>(strs[i].size()));
            }
        }
        return ans;
    }
} // namespace findLUSlength2

void findLUSlength2_test() {
    vector<string> strs = {"aba", "cdc", "eae"};
    cout << "字符串数组：";
    print_vector(strs);
    cout << " 的最长特殊序列个数为 " << findLUSlength2::findLUSlength(strs)
         << endl;
    strs = {"aaa", "aaa", "aa"};
    cout << "字符串数组：";
    print_vector(strs);
    cout << " 的最长特殊序列个数为 " << findLUSlength2::findLUSlength(strs)
         << endl;
}

namespace checkSubarraySum {
    bool checkSubarraySum(vector<int> &nums, int k) {
        int m = nums.size();
        if (m < 2) {
            return false;
        }
        unordered_map<int, int> mp;
        mp[0] = -1;
        int remainder = 0;
        for (int i = 0; i < m; i++) {
            remainder = (remainder + nums[i]) % k;
            if (mp.count(remainder)) {
                int prevIndex = mp[remainder];
                if (i - prevIndex >= 2) {
                    return true;
                }
            } else {
                mp[remainder] = i;
            }
        }
        return false;
    }
} // namespace checkSubarraySum

void checkSubarraySum_test() {
    vector<int> nums;
    int k;
    nums = {23, 2, 4, 6, 7};
    k = 6;
    cout << "数组 ";
    print_vector(nums);
    if (checkSubarraySum::checkSubarraySum(nums, k)) {
        cout << "存在和为 " << k << " 的连续子数组" << endl;
    } else {
        cout << "不存在和为 " << k << " 的连续子数组" << endl;
    }
    k = 6;
    cout << "数组 ";
    print_vector(nums);
    if (checkSubarraySum::checkSubarraySum(nums, k)) {
        cout << "存在和为 " << k << " 的连续子数组" << endl;
    } else {
        cout << "不存在和为 " << k << " 的连续子数组" << endl;
    }
    k = 13;
    cout << "数组 ";
    print_vector(nums);
    if (checkSubarraySum::checkSubarraySum(nums, k)) {
        cout << "存在和为 " << k << " 的连续子数组" << endl;
    } else {
        cout << "不存在和为 " << k << " 的连续子数组" << endl;
    }
}

namespace findMaxLength {
    int findMaxLength(vector<int> &nums) {
        int maxLength = 0;
        unordered_map<int, int> mp;
        int counter = 0;
        mp[counter] = -1;
        int n = nums.size();
        for (int i = 0; i < n; i++) {
            int num = nums[i];
            if (num == 1) {
                counter++;
            } else {
                counter--;
            }
            if (mp.count(counter)) {
                int prevIndex = mp[counter];
                maxLength = max(maxLength, i - prevIndex);
            } else {
                mp[counter] = i;
            }
        }
        return maxLength;
    }
}

void findMaxLength_test() {
    vector<int> nums;
    nums = {0, 1};
    cout << findMaxLength::findMaxLength(nums) << endl;
    nums = {0, 1, 0};
    cout << findMaxLength::findMaxLength(nums) << endl;
}

namespace leastBricks {
    int leastBricks(vector<vector<int>> &wall) {
        unordered_map<int, int> cnt;
        for (auto &widths : wall) {
            int n = widths.size();
            int sum = 0;
            for (int i = 0; i < n - 1; ++i) {
                sum += widths[i];
                cnt[sum]++;
            }
        }
        int maxCnt = 0;
        for (auto[_, c] : cnt) {
            maxCnt = max(maxCnt, c);
        }
        return wall.size() - maxCnt;
    }
}

void leastBricks_test() {
    vector<vector<int>> wall;
    wall = {{1, 2, 2, 1},
            {3, 1, 2},
            {1, 3, 2},
            {2, 4},
            {3, 1, 2},
            {1, 3, 1, 1}};
    cout << leastBricks::leastBricks(wall) << endl;
    wall = {{1},
            {1},
            {1}};
    cout << leastBricks::leastBricks(wall) << endl;
}

namespace nextGreaterElement3 {
    int nextGreaterElement(int n) {
        string nums = to_string(n);
        int i = nums.size() - 2;
        while (i >= 0 && nums[i] >= nums[i + 1]) {
            i--;
        }
        if (i < 0) {
            return -1;
        }
        int j = nums.size() - 1;
        while (j >= 0 && nums[i] >= nums[j]) {
            j--;
        }
        swap(nums[i], nums[j]);
        reverse(nums.begin() + i + 1, nums.end());
        long ans = stol(nums);
        return ans > INT_MAX ? -1 : ans;
    }
}

void nextGreaterElement3_test() {
    int n;
    n = 12;
    cout << nextGreaterElement3::nextGreaterElement(n) << endl;
    n = 21;
    cout << nextGreaterElement3::nextGreaterElement(n) << endl;
}

namespace QTree {
    class Node {
    public:
        bool val;
        bool isLeaf;
        Node *topLeft;
        Node *topRight;
        Node *bottomLeft;
        Node *bottomRight;

        Node() {
            val = false;
            isLeaf = false;
            topLeft = NULL;
            topRight = NULL;
            bottomLeft = NULL;
            bottomRight = NULL;
        }

        Node(bool _val, bool _isLeaf) {
            val = _val;
            isLeaf = _isLeaf;
            topLeft = NULL;
            topRight = NULL;
            bottomLeft = NULL;
            bottomRight = NULL;
        }

        Node(bool _val, bool _isLeaf, Node *_topLeft, Node *_topRight, Node *_bottomLeft, Node *_bottomRight) {
            val = _val;
            isLeaf = _isLeaf;
            topLeft = _topLeft;
            topRight = _topRight;
            bottomLeft = _bottomLeft;
            bottomRight = _bottomRight;
        }

    };

    Node *create_qtreenode(vector<vector<int>> &qtree) {
        Node *node = new Node(qtree[0][1], qtree[0][0]);
        queue<Node *> q;
        q.push(node);
        int start = 0;
        while (!q.empty()) {
            int size = q.size();
            start = start + size;
            for (int i = 0; i < size; ++i) {
                auto tmp = q.front();
                if (tmp->isLeaf == 0) {
                    tmp->topLeft = new Node(qtree[start + 0][1], qtree[start + 0][0]);
                    tmp->topRight = new Node(qtree[start + 1][1], qtree[start + 1][0]);
                    tmp->bottomLeft = new Node(qtree[start + 2][1], qtree[start + 2][0]);
                    tmp->bottomRight = new Node(qtree[start + 3][1], qtree[start + 3][0]);
                    q.push(tmp->topLeft);
                    q.push(tmp->topRight);
                    q.push(tmp->bottomLeft);
                    q.push(tmp->bottomRight);
                } else if (tmp->isLeaf == 1) {
                    start = start + 4;
                }
                q.pop();
            }
        }
        return node;
    }

    Node *intersect(Node *quadTree1, Node *quadTree2) {
        if (quadTree1->isLeaf) {
            if (quadTree1->val) {
                return new Node(true, true);
            }
            return new Node(quadTree2->val, quadTree2->isLeaf, quadTree2->topLeft, quadTree2->topRight,
                            quadTree2->bottomLeft, quadTree2->bottomRight);
        }
        if (quadTree2->isLeaf) {
            return intersect(quadTree2, quadTree1);
        }
        Node *o1 = intersect(quadTree1->topLeft, quadTree2->topLeft);
        Node *o2 = intersect(quadTree1->topRight, quadTree2->topRight);
        Node *o3 = intersect(quadTree1->bottomLeft, quadTree2->bottomLeft);
        Node *o4 = intersect(quadTree1->bottomRight, quadTree2->bottomRight);
        if (o1->isLeaf && o2->isLeaf && o3->isLeaf && o4->isLeaf && o1->val == o2->val && o1->val == o3->val &&
            o1->val == o4->val) {
            return new Node(o1->val, true);
        }
        return new Node(false, false, o1, o2, o3, o4);
    }

    void print_qtreenode(Node *root) {

    }
}

void QTree_test() {
    vector<vector<int>> qtree1, qtree2;
    qtree1 = {{0, 1},
              {1, 1},
              {1, 1},
              {1, 0},
              {1, 0}};
    qtree2 = {{0,  1},
              {1,  1},
              {0,  1},
              {1,  1},
              {1,  0},
              {-1, -1},
              {-1, -1},
              {-1, -1},
              {-1, -1},
              {1,  0},
              {1,  0},
              {1,  1},
              {1,  1}};
    QTree::Node *qtree_node1 = QTree::create_qtreenode(qtree1);
    QTree::Node *qtree_node2 = QTree::create_qtreenode(qtree2);
    QTree::Node *res = QTree::intersect(qtree_node1, qtree_node2);
    qtree1 = {{1, 0}};
    qtree2 = {{1, 0}};
    qtree_node1 = QTree::create_qtreenode(qtree1);
    qtree_node2 = QTree::create_qtreenode(qtree2);
    res = QTree::intersect(qtree_node1, qtree_node2);
}

namespace NTree {
    class Node {
    public:
        int val;
        vector<Node *> children;

        Node() {}

        Node(int _val) {
            val = _val;
        }

        Node(int _val, vector<Node *> _children) {
            val = _val;
            children = _children;
        }
    };

    Node *createNTree(vector<int> tree) {
        Node *root = new Node(tree[0]);
        queue<Node *> q;
        int start = 2;
        q.push(root);
        while (!q.empty()) {
            int size = q.size();

            for (int i = 0; i < size; ++i) {
                auto n = q.front();
                while (start < tree.size() && tree[start] != 0) {
                    auto tmp = new Node(tree[start]);
                    n->children.push_back(tmp);
                    q.push(tmp);
                    start++;
                }
                if (start < tree.size() && tree[start] == 0) {
                    start++;
                }
                q.pop();
            }
        }
        return root;
    }

    int maxDepth(Node *root) {
        queue<Node *> q;
        q.push(root);
        int res = 0;
        if (root == nullptr) {
            return 0;
        }
        int laycount = 0;
        int laynum = 1;
        while (!q.empty()) {
            auto tmp = q.front();
            for (auto node : tmp->children) {
                q.push(node);
            }
            q.pop();
            laycount++;
            if (laycount == laynum) {
                res++;
                laynum = q.size();
                laycount = 0;
            }
        }
        return res;
    }
}

void ntreedepth_test() {
    vector<int> ntree = {1, 0, 3, 2, 4, 0, 5, 6};
    NTree::Node *node = NTree::createNTree(ntree);
    cout << NTree::maxDepth(node) << endl;
    ntree = {1, 0, 2, 3, 4, 5, 0, 0, 6, 7, 0, 8, 0, 9, 10, 0, 0, 11, 0, 12, 0, 13, 0, 0, 14};
    node = NTree::createNTree(ntree);
    cout << NTree::maxDepth(node) << endl;
}

namespace subarraySum {
    int subarraySum(vector<int> &nums, int k) {
        unordered_map<int, int> mp;
        mp[0] = 1;
        int count = 0, pre = 0;
        for (auto x : nums) {
            pre += x;
            if (mp.find(pre - k) != mp.end()) {
                count += mp[pre - k];
            }
            mp[pre]++;
        }
        return count;
    }
}

void subarraySum_test() {
    vector<int> nums;
    int k;
    nums = {1, 1, 1};
    k = 2;
    cout << subarraySum::subarraySum(nums, k) << endl;
    nums = {1, 2, 3};
    k = 3;
    cout << subarraySum::subarraySum(nums, k) << endl;
}

namespace arrayPairSum {
    int arrayPairSum(vector<int> &nums) {
        sort(nums.begin(), nums.end());
        int ans = 0;
        for (int i = 0; i < nums.size(); i += 2) {
            ans += nums[i];
        }
        return ans;
    }
}

void arrayPairSum_test() {
    vector<int> nums;
    nums = {1, 4, 3, 2};
    cout << arrayPairSum::arrayPairSum(nums) << endl;
    nums = {6, 2, 6, 5, 1, 2};
    cout << arrayPairSum::arrayPairSum(nums) << endl;
}

namespace findTilt {
    int ans = 0;

    int dfs(TreeNode::TreeNode *node) {
        if (node == nullptr) {
            return 0;
        }
        int left = dfs(node->left);
        int righ = dfs(node->right);
        ans += abs(left - righ);
        return left + righ + node->val;
    }

    int findTilt(TreeNode::TreeNode *root) {
        ans = 0;
        dfs(root);
        return ans;
    }
}

void findTilt_test() {
    TreeNode::TreeNode *node;
    vector<int> tree;
    tree = {1, 2, 3};
    node = create_treenode(tree);
    cout << findTilt::findTilt(node) << endl;
    tree = {4, 2, 9, 3, 5, 0, 7};
    node = create_treenode(tree);
    cout << findTilt::findTilt(node) << endl;
    tree = {21, 7, 14, 1, 1, 2, 2, 3, 3};
    node = create_treenode(tree);
    cout << findTilt::findTilt(node) << endl;
}

namespace nearestPalindromic {
    using ULL = unsigned long long;

    vector<ULL> getCandidates(const string &n) {
        int len = n.length();
        vector<ULL> candidates = {
                (ULL) pow(10, len - 1) - 1,
                (ULL) pow(10, len) + 1,
        };
        ULL selfPrefix = stoull(n.substr(0, (len + 1) / 2));
        for (int i : {selfPrefix - 1, selfPrefix, selfPrefix + 1}) {
            string prefix = to_string(i);
            string candidate = prefix + string(prefix.rbegin() + (len & 1), prefix.rend());
            candidates.push_back(stoull(candidate));
        }
        return candidates;
    }

    string nearestPalindromic(string n) {
        ULL selfNumber = stoull(n), ans = -1;
        const vector<ULL> &candidates = getCandidates(n);
        for (auto &candidate : candidates) {
            if (candidate != selfNumber) {
                if (ans == -1 ||
                    llabs(candidate - selfNumber) < llabs(ans - selfNumber) ||
                    llabs(candidate - selfNumber) == llabs(ans - selfNumber) && candidate < ans) {
                    ans = candidate;
                }
            }
        }
        return to_string(ans);
    }
}

void nearestPalindromic_test() {
    string n;
    n = "123";
    cout << n << " 的最近的回文数是 " << nearestPalindromic::nearestPalindromic(n) << endl;
    n = "1";
    cout << n << " 的最近的回文数是 " << nearestPalindromic::nearestPalindromic(n) << endl;
    n = "99321";
    cout << n << " 的最近的回文数是 " << nearestPalindromic::nearestPalindromic(n) << endl;
}

namespace arrayNesting {
    int arrayNesting(vector<int> &nums) {
        int ans = 0, n = nums.size();
        vector<int> vis(n);
        for (int i = 0; i < n; ++i) {
            int cnt = 0;
            while (!vis[i]) {
                vis[i] = true;
                i = nums[i];
                ++cnt;
            }
            ans = max(ans, cnt);
        }
        return ans;
    }
}

void arrayNesting_test() {
    vector<int> n;
    n = {5, 4, 0, 3, 1, 6, 2};
    cout << arrayNesting::arrayNesting(n) << endl;
    n = {0, 1, 2};
    cout << arrayNesting::arrayNesting(n) << endl;
}

namespace matrixReshape {
    vector<vector<int>> matrixReshape(vector<vector<int>> &mat, int r, int c) {
        int m = mat.size();
        int n = mat[0].size();
        if (m * n != r * c) {
            return mat;
        }

        vector<vector<int>> ans(r, vector<int>(c));
        for (int x = 0; x < m * n; ++x) {
            ans[x / c][x % c] = mat[x / n][x % n];
        }
        return ans;

    }
}

void matrixReshape_test() {
    vector<vector<int>> mat;
    int r, c;
    mat = {{1, 2},
           {3, 4}};
    r = 1, c = 4;
    auto ans = matrixReshape::matrixReshape(mat, r, c);
    print_mat(ans);
    r = 2, c = 4;
    ans = matrixReshape::matrixReshape(mat, r, c);
    print_mat(ans);
}

namespace checkInclusion {
    bool checkInclusion(string s1, string s2) {
        int n = s1.length(), m = s2.length();
        if (n > m) {
            return false;
        }
        vector<int> cnt1(26), cnt2(26);
        for (int i = 0; i < n; ++i) {
            ++cnt1[s1[i] - 'a'];
            ++cnt2[s2[i] - 'a'];
        }
        if (cnt1 == cnt2) {
            return true;
        }
        for (int i = n; i < m; ++i) {
            ++cnt2[s2[i] - 'a'];
            --cnt2[s2[i - n] - 'a'];
            if (cnt1 == cnt2) {
                return true;
            }
        }
        return false;

    }
}

void checkInclusion_test() {
    string s1, s2;
    s1 = "ab";
    s2 = "eidbaooo";
    cout << checkInclusion::checkInclusion(s1, s2) << endl;
    s1 = "ab";
    s2 = "eidboaoo";
    cout << checkInclusion::checkInclusion(s1, s2) << endl;
}

namespace isSubtree {
    bool check(TreeNode::TreeNode *o, TreeNode::TreeNode *t) {
        if (!o && !t) {
            return true;
        }
        if ((o && !t) || (!o && t) || (o->val != t->val)) {
            return false;
        }
        return check(o->left, t->left) && check(o->right, t->right);
    }

    bool dfs(TreeNode::TreeNode *o, TreeNode::TreeNode *t) {
        if (!o) {
            return false;
        }
        return check(o, t) || dfs(o->left, t) || dfs(o->right, t);
    }

    bool isSubtree(TreeNode::TreeNode *root, TreeNode::TreeNode *subRoot) {
        return dfs(root, subRoot);
    }
}

void isSubtree_test() {
    TreeNode::TreeNode *root, *subRoot;
    vector<int> root_data, subroot_data;
    root_data = {3, 4, 5, 1, 2};
    subroot_data = {4, 1, 2};
    root = create_treenode(root_data);
    subRoot = create_treenode(subroot_data);
    cout << isSubtree::isSubtree(root, subRoot);
    root_data = {3, 4, 5, 1, 2, 0, 0, 0, 0, 10};
    subroot_data = {4, 1, 2};
    cout << isSubtree::isSubtree(root, subRoot);
}

namespace distributeCandies {
    int distributeCandies(vector<int> &candyType) {
        unordered_set<int> set(candyType.begin(), candyType.end());

        return min(set.size(), candyType.size() / 2);
    }
}

void distributeCandies_test() {
    vector<int> candyType;
    candyType = {1, 1, 2, 2, 3, 3};
    cout << distributeCandies::distributeCandies(candyType) << endl;
    candyType = {1, 1, 2, 3};
    cout << distributeCandies::distributeCandies(candyType) << endl;
    candyType = {6, 6, 6, 6};
    cout << distributeCandies::distributeCandies(candyType) << endl;
}

namespace findPaths {
    static constexpr int MOD = 1'000'000'007;

    int findPaths(int m, int n, int maxMove, int startRow, int startColumn) {
        vector<vector<int>> directions = {{-1, 0},
                                          {1,  0},
                                          {0,  -1},
                                          {0,  1}};
        int outCounts = 0;
        vector<vector<int>> dp(m, vector<int>(n));
        dp[startRow][startColumn] = 1;
        for (int i = 0; i < maxMove; ++i) {
            vector<vector<int>> dpNew(m, vector<int>(n));
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < n; ++k) {
                    int count = dp[j][k];
                    if (count > 0) {
                        for (auto direction: directions) {
                            int j1 = j + direction[0];
                            int k1 = k + direction[1];
                            if (j1 >= 0 && j1 < m && k1 >= 0 && k1 < n) {
                                dpNew[j1][k1] = (dpNew[j1][k1] + count) % MOD;
                            } else {
                                outCounts = (outCounts + count) % MOD;
                            }
                        }
                    }
                }
            }
            dp = dpNew;
        }
        return outCounts;
    }
}

void findPaths_test() {
    int m = 2, n = 2, maxMove = 2, startRow = 0, startColumn = 0;
    cout << findPaths::findPaths(m, n, maxMove, startRow, startColumn) << endl;
    m = 1, n = 3, maxMove = 3, startRow = 0, startColumn = 1;
    cout << findPaths::findPaths(m, n, maxMove, startRow, startColumn) << endl;
}

namespace findUnsortedSubarray {
    int findUnsortedSubarray(vector<int> &nums) {
        int srt = 0;
        int end = nums.size() - 1;
        auto sort_nums = nums;
        sort(sort_nums.begin(), sort_nums.end());

        while (true) {
            if (srt >= end) {
                return 0;

            }
            if (sort_nums[srt] != nums[srt] && sort_nums[end] != nums[end])
                break;
            if (sort_nums[srt] == nums[srt])
                srt++;
            if (sort_nums[end] == nums[end])
                end--;
        }
        return end - srt + 1;
    }
}

void findUnsortedSubarray_test() {
    vector<int> nums;
    nums = {2, 6, 4, 8, 10, 9, 15};
    cout << findUnsortedSubarray::findUnsortedSubarray(nums) << endl;
    nums = {1, 2, 3, 4};
    cout << findUnsortedSubarray::findUnsortedSubarray(nums) << endl;
}

namespace minDistance {
    int minDistance(string word1, string word2) {
        int m = word1.size();
        int n = word2.size();
        vector<vector<int>> dp(m + 1, vector<int>(n + 1));

        for (int i = 1; i <= m; i++) {
            char c1 = word1[i - 1];
            for (int j = 1; j <= n; j++) {
                char c2 = word2[j - 1];
                if (c1 == c2) {
                    dp[i][j] = dp[i - 1][j - 1] + 1;
                } else {
                    dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
                }
            }
        }

        int lcs = dp[m][n];
        return m - lcs + n - lcs;
    }
}

void minDistance_test() {
    string word1, word2;
    word1 = "sea";
    word2 = "eat";
    cout << minDistance::minDistance(word1, word2) << endl;
    word1 = "leetcode";
    word2 = "etco";
    cout << minDistance::minDistance(word1, word2) << endl;
}

namespace outerTrees {
    int cross(vector<int> &p, vector<int> &q, vector<int> &r) {
        return (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]);
    }

    vector<vector<int>> outerTrees(vector<vector<int>> &trees) {
        int n = trees.size();
        if (n < 4) {
            return trees;
        }
        int leftMost = 0;
        for (int i = 0; i < n; i++) {
            if (trees[i][0] < trees[leftMost][0] ||
                (trees[i][0] == trees[leftMost][0] &&
                 trees[i][1] < trees[leftMost][1])) {
                leftMost = i;
            }
        }

        vector<vector<int>> res;
        vector<bool> visit(n, false);
        int p = leftMost;
        do {
            int q = (p + 1) % n;
            for (int r = 0; r < n; r++) {
                /* 如果 r 在 pq 的右侧，则 q = r */
                if (cross(trees[p], trees[q], trees[r]) < 0) {
                    q = r;
                }
            }
            /* 是否存在点 i, 使得 p 、q 、i 在同一条直线上 */
            for (int i = 0; i < n; i++) {
                if (visit[i] || i == p || i == q) {
                    continue;
                }
                if (cross(trees[p], trees[q], trees[i]) == 0) {
                    res.emplace_back(trees[i]);
                    visit[i] = true;
                }
            }
            if (!visit[q]) {
                res.emplace_back(trees[q]);
                visit[q] = true;
            }
            p = q;
        } while (p != leftMost);
        return res;
    }
}

void outerTrees_test() {
    vector<vector<int>> trees = {{1, 1},
                                 {2, 2},
                                 {2, 0},
                                 {2, 4},
                                 {3, 3},
                                 {4, 2}};
    auto ans = outerTrees::outerTrees(trees);
    for (auto line : ans) {
        print_vector(line);
    }
    trees = {{1, 2},
             {2, 2},
             {4, 2}};
    ans = outerTrees::outerTrees(trees);
    for (auto line : ans) {
        print_vector(line);
    }
}

namespace preorder {
    void dfs(NTree::Node *root, vector<int> &ans) {
        if (root == nullptr)
            return;
        ans.push_back(root->val);
        for (auto child : root->children) {
            dfs(child, ans);
        }
    }

    vector<int> preorder(NTree::Node *root) {
        vector<int> ans;
        dfs(root, ans);
        return ans;
    }
}

void preorder_test() {
    vector<int> ans;
    vector<int> tree = {1, 0, 3, 2, 4, 0, 5, 6};
    NTree::Node *root = NTree::createNTree(tree);
    ans = preorder::preorder(root);
    print_vector(ans);
    tree = {1, 0, 2, 3, 4, 5, 0, 0, 6, 7, 0, 8, 0, 9, 10, 0, 0, 11, 0, 12, 0, 13, 0, 0, 14};
    root = NTree::createNTree(tree);
    ans = preorder::preorder(root);
    print_vector(ans);
}

namespace postorder {
    void dfs(NTree::Node *root, vector<int> &ans) {
        if (root == nullptr)
            return;

        for (auto child : root->children) {
            dfs(child, ans);
        }
        ans.push_back(root->val);
    }

    vector<int> postorder(NTree::Node *root) {
        vector<int> ans;
        dfs(root, ans);
        return ans;
    }
}

void postorder_test() {
    vector<int> ans;
    vector<int> tree = {1, 0, 3, 2, 4, 0, 5, 6};
    NTree::Node *root = NTree::createNTree(tree);
    ans = postorder::postorder(root);
    print_vector(ans);
    tree = {1, 0, 2, 3, 4, 5, 0, 0, 6, 7, 0, 8, 0, 9, 10, 0, 0, 11, 0, 12, 0, 13, 0, 0, 14};
    root = NTree::createNTree(tree);
    ans = postorder::postorder(root);
    print_vector(ans);
}

namespace fractionAddition {
    string fractionAddition(string expression) {
        long long x = 0, y = 1; // 分子，分母
        int index = 0, n = expression.size();
        while (index < n) {
            // 读取分子
            long long x1 = 0, sign = 1;
            if (expression[index] == '-' || expression[index] == '+') {
                sign = expression[index] == '-' ? -1 : 1;
                index++;
            }
            while (index < n && isdigit(expression[index])) {
                x1 = x1 * 10 + expression[index] - '0';
                index++;
            }
            x1 = sign * x1;
            index++;

            // 读取分母
            long long y1 = 0;
            while (index < n && isdigit(expression[index])) {
                y1 = y1 * 10 + expression[index] - '0';
                index++;
            }

            x = x * y1 + x1 * y;
            y *= y1;
        }
        if (x == 0) {
            return "0/1";
        }
        long long g = gcd(abs(x), y); // 获取最大公约数
        return to_string(x / g) + "/" + to_string(y / g);
    }
}

void fractionAddition_test() {
    string expression = "-1/2+1/2";
    cout << expression << " = " << fractionAddition::fractionAddition(expression) << endl;
    expression = "-1/2+1/2+1/3";
    cout << expression << " = " << fractionAddition::fractionAddition(expression) << endl;
    expression = "1/3-1/2";
    cout << expression << " = " << fractionAddition::fractionAddition(expression) << endl;
}

namespace findLHS {
    int findLHS(vector<int> &nums) {
        unordered_map<int, int> cnt;
        int res = 0;
        for (int num : nums) {
            cnt[num]++;
        }
        for (auto[key, val] : cnt) {
            if (cnt.count(key + 1)) {
                res = max(res, val + cnt[key + 1]);
            }
        }
        return res;
    }
}

void findLHS_test() {
    vector<int> nums;
    nums = {1, 3, 2, 2, 5, 2, 3, 7};
    cout << findLHS::findLHS(nums) << endl;
    nums = {1, 2, 3, 4};
    cout << findLHS::findLHS(nums) << endl;
    nums = {1, 1, 1, 1};
    cout << findLHS::findLHS(nums) << endl;
}

namespace triangleNumber {
    int triangleNumber(vector<int> &nums) {
        int n = nums.size();
        sort(nums.begin(), nums.end());
        int ans = 0;
        for (int i = 0; i < n; ++i) {
            int k = i;
            for (int j = i + 1; j < n; ++j) {
                while (k + 1 < n && nums[k + 1] < nums[i] + nums[j]) {
                    ++k;
                }
                ans += max(k - j, 0);
            }
        }
        return ans;
    }
}

void triangleNumber_test() {
    vector<int> nums;
    nums = {2, 2, 3, 4};
    cout << triangleNumber::triangleNumber(nums) << endl;
    nums = {4, 2, 3, 4};
    cout << triangleNumber::triangleNumber(nums) << endl;
    nums = {48, 66, 61, 46, 94, 75};
    cout << triangleNumber::triangleNumber(nums) << endl;
}

namespace mergeTrees {
    TreeNode::TreeNode *mergeTrees(TreeNode::TreeNode *root1, TreeNode::TreeNode *root2) {
        if (root1 == nullptr) {
            return root2;
        }
        if (root2 == nullptr) {
            return root1;
        }
        auto merged = new TreeNode::TreeNode(root1->val + root2->val);
        merged->left = mergeTrees(root1->left, root2->left);
        merged->right = mergeTrees(root1->right, root2->right);
        return merged;

    }
}

void mergeTrees_test() {
    vector<int> nums1, nums2;
    nums1 = {1, 3, 2, 5};
    nums2 = {2, 1, 3, 0, 4, 0, 7};
    TreeNode::TreeNode *root1 = create_treenode(nums1);
    TreeNode::TreeNode *root2 = create_treenode(nums2);
    auto ans = mergeTrees::mergeTrees(root1, root2);
    cout << TreeNode::print_tree(ans) << endl;
    nums1 = {1};
    nums2 = {1, 2};
    root1 = create_treenode(nums1);
    root2 = create_treenode(nums2);
    ans = mergeTrees::mergeTrees(root1, root2);
    cout << TreeNode::print_tree(ans) << endl;
}

namespace leastInterval {
    int leastInterval(vector<char> &tasks, int n) {
        unordered_map<char, int> freq;
        for (char ch: tasks) {
            ++freq[ch];
        }

        // 最多的执行次数
        int maxExec = max_element(freq.begin(), freq.end(), [](const auto &u, const auto &v) {
            return u.second < v.second;
        })->second;
        // 具有最多执行次数的任务数量
        int maxCount = accumulate(freq.begin(), freq.end(), 0, [=](int acc, const auto &u) {
            return acc + (u.second == maxExec);
        });

        return max((maxExec - 1) * (n + 1) + maxCount, static_cast<int>(tasks.size()));
    }
}

void leastInterval_test() {
    vector<char> tasks;
    int n = 2;
    tasks = {'A', 'A', 'A', 'B', 'B', 'B'};// 8
    cout << leastInterval::leastInterval(tasks, n) << endl;
    n = 1;
    tasks = {'A', 'C', 'A', 'B', 'D', 'B'};//6
    cout << leastInterval::leastInterval(tasks, n) << endl;
    n = 3;
    tasks = {'A', 'A', 'A', 'B', 'B', 'B'};//10
    cout << leastInterval::leastInterval(tasks, n) << endl;
    n = 0;
    tasks = {'A', 'A', 'A', 'B', 'B', 'B'};//6
    cout << leastInterval::leastInterval(tasks, n) << endl;
    n = 2;
    tasks = {'A', 'A', 'A', 'A', 'A', 'A', 'B', 'C', 'D', 'E', 'F', 'G'};//16
    cout << leastInterval::leastInterval(tasks, n) << endl;
}

namespace addOneRow {
    TreeNode::TreeNode *addOneRow(TreeNode::TreeNode *root, int val, int depth) {
        if (root == nullptr) {
            return nullptr;
        }
        if (depth == 1) {
            return new TreeNode::TreeNode(val, root, nullptr);
        }
        if (depth == 2) {
            root->left = new TreeNode::TreeNode(val, root->left, nullptr);
            root->right = new TreeNode::TreeNode(val, nullptr, root->right);
        } else {
            root->left = addOneRow(root->left, val, depth - 1);
            root->right = addOneRow(root->right, val, depth - 1);
        }
        return root;
    }
}

void addOneRow_test() {
    vector<int> tree;
    int val, depth;
    val = 1;
    depth = 2;
    tree = {4, 2, 6, 3, 1, 5};
    auto root = create_treenode(tree);
    auto ans = addOneRow::addOneRow(root, val, depth);
    cout << TreeNode::print_tree(ans) << endl;
    val = 1;
    depth = 3;
    tree = {4, 2, 0, 3, 1};
    root = create_treenode(tree);
    ans = addOneRow::addOneRow(root, val, depth);
    cout << TreeNode::print_tree(ans) << endl;
}

namespace maximumProduct {
    int maximumProduct(vector<int> &nums) {
        sort(nums.begin(), nums.end());
        int ans = 0;
        int a, b, c;
        a = nums[nums.size() - 3];
        b = nums[nums.size() - 2];
        c = nums[nums.size() - 1];
        int d, e, f;
        d = nums[0];
        e = nums[1];
        f = nums[2];
        if (d * e * c > a * b * c) {
            a = d;
            b = e;
        }
        return a * b * c;
    }
}

void maximumProduct_test() {
    vector<int> nums;
    nums = {1, 2, 3};
    cout << maximumProduct::maximumProduct(nums) << endl;
    nums = {1, 2, 3, 4};
    cout << maximumProduct::maximumProduct(nums) << endl;
    nums = {-1, -2, -3};
    cout << maximumProduct::maximumProduct(nums) << endl;
    nums = {-100, -98, -1, 2, 3, 4};
    cout << maximumProduct::maximumProduct(nums) << endl;
    nums = {-8, -7, -2, 10, 20};
    cout << maximumProduct::maximumProduct(nums) << endl;
}

namespace kInversePairs {
    static constexpr int mod = 1000000007;

    int kInversePairs(int n, int k) {
        vector<vector<int>> f(2, vector<int>(k + 1));
        f[0][0] = 1;
        for (int i = 1; i <= n; ++i) {
            for (int j = 0; j <= k; ++j) {
                int cur = i & 1, prev = cur ^1;
                f[cur][j] = (j - 1 >= 0 ? f[cur][j - 1] : 0) - (j - i >= 0 ? f[prev][j - i] : 0) + f[prev][j];
                if (f[cur][j] >= mod) {
                    f[cur][j] -= mod;
                } else if (f[cur][j] < 0) {
                    f[cur][j] += mod;
                }
            }
        }
        return f[n & 1][k];
    }
}

void kInversePairs_test() {
    int n, k;
    n = 3, k = 0;
    cout << kInversePairs::kInversePairs(n, k) << endl;
    n = 3, k = 1;
    cout << kInversePairs::kInversePairs(n, k) << endl;
}

namespace judgeSquareSum {
    bool judgeSquareSum(int c) {
        long left = 0;
        long right = (long) sqrt(c);
        while (left <= right) {
            long sum = left * left + right * right;
            if (sum == c) {
                return true;
            } else if (sum > c) {
                right--;
            } else {
                left++;
            }
        }
        return false;

    }
}

void judgeSquareSum_test() {
    int c;
    c = 5;
    cout << "c:" << c << " 是平方数之和，" << judgeSquareSum::judgeSquareSum(c) << endl;
    c = 3;
    cout << "c:" << c << " 是平方数之和，" << judgeSquareSum::judgeSquareSum(c) << endl;
    c = 4;
    cout << "c:" << c << " 是平方数之和，" << judgeSquareSum::judgeSquareSum(c) << endl;
}

namespace averageOfLevels {
    vector<double> averageOfLevels(TreeNode::TreeNode *root) {
        auto averages = vector<double>();
        auto q = queue<TreeNode::TreeNode *>();
        q.push(root);
        while (!q.empty()) {
            double sum = 0;
            int size = q.size();
            for (int i = 0; i < size; i++) {
                auto node = q.front();
                q.pop();
                sum += node->val;
                auto left = node->left, right = node->right;
                if (left != nullptr) {
                    q.push(left);
                }
                if (right != nullptr) {
                    q.push(right);
                }
            }
            averages.push_back(sum / size);
        }
        return averages;
    }
}

void averageOfLevels_test() {
    vector<int> nums;
    nums = {3, 9, 20, 0, 0, 15, 7};
    TreeNode::TreeNode *root = create_treenode(nums);
    auto ans = averageOfLevels::averageOfLevels(root);
    print_vector(ans);
    nums = {3, 9, 20, 15, 7};
    root = create_treenode(nums);
    ans = averageOfLevels::averageOfLevels(root);
    print_vector(ans);
}

namespace findErrorNums {
    vector<int> findErrorNums(vector<int> &nums) {
        vector<int> ans;
        unordered_set<int> set;
        int loss = 0;
        for (int i = 0; i < nums.size(); ++i) {
            if (set.find(nums[i]) != set.end()) {
                ans.push_back(nums[i]);
            } else {
                set.insert(nums[i]);
            }
        }
        for (int i = 1; i <= nums.size(); ++i) {
            if (set.find(i) == set.end()) {
                ans.push_back(i);
            }
        }
        return ans;
    }
}

void findErrorNums_test() {
    vector<int> nums;
    nums = {1, 2, 2, 4};
    auto ans = findErrorNums::findErrorNums(nums);
    print_vector(ans);
    nums = {1, 1};
    ans = findErrorNums::findErrorNums(nums);
    print_vector(ans);
    nums = {2, 2};
    ans = findErrorNums::findErrorNums(nums);
    print_vector(ans);
}

namespace findLongestChain {
    int findLongestChain(vector<vector<int>> &pairs) {
        int n = pairs.size();
        sort(pairs.begin(), pairs.end());
        vector<int> dp(n, 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (pairs[i][0] > pairs[j][1]) {
                    dp[i] = max(dp[i], dp[j] + 1);
                }
            }
        }
        return dp[n - 1];
    }
}

void findLongestChain_test() {
    vector<vector<int>> pairs;
    pairs = {{1, 2},
             {2, 3},
             {3, 4}};
    cout << findLongestChain::findLongestChain(pairs) << endl;
    pairs = {{1, 2},
             {7, 8},
             {4, 5}};
    cout << findLongestChain::findLongestChain(pairs) << endl;
    pairs = {{1, 2}};
    cout << findLongestChain::findLongestChain(pairs) << endl;
}

namespace countSubstrings {
    int countSubstrings(string s) {
        int n = s.size(), ans = 0;
        for (int i = 0; i < 2 * n - 1; ++i) {
            int l = i / 2, r = i / 2 + i % 2;
            while (l >= 0 && r < n && s[l] == s[r]) {
                --l;
                ++r;
                ++ans;
            }
        }
        return ans;
    }
}

void countSubstrings_test() {
    string s;
    s = "abc";
    cout << s << "的回文子串数：" << countSubstrings::countSubstrings(s) << std::endl;
    s = "aaa";
    cout << s << "的回文子串数：" << countSubstrings::countSubstrings(s) << std::endl;
}


#include <regex>

namespace replaceWords {
    std::vector<std::string> split(const std::string &str, const std::string &delimiter) {
        std::regex re(delimiter);
        std::sregex_token_iterator it(str.begin(), str.end(), re, -1);
        std::sregex_token_iterator reg_end;

        std::vector<std::string> result(it, reg_end);
        return result;
    }

    string replaceWords(vector<string> &dictionary, string sentence) {
        string ans;
        unordered_set<string> set;
        for (auto s:dictionary) {
            set.insert(s);
        }
        vector<string> tokens;
        tokens = split(sentence, "\\s+");
        for (int j = 0; j < tokens.size(); ++j) {
            auto t = tokens[j];
            bool flag = false;
            for (int i = 0; i < t.size(); ++i) {
                string tmp;
                tmp.assign(t.begin(), t.begin() + i);
                if (set.find(tmp) != set.end()) {
                    ans += tmp;
                    flag = true;
                    break;
                }
            }
            if (!flag) {
                ans += t;
            }
            if (j < tokens.size() - 1)
                ans += " ";

        }
        return ans;
    }
}

void replaceWords_test() {
    vector<string> dictionary;
    string sentence;
    dictionary = {"cat", "bat", "rat"};
    sentence = "the cattle was rattled by the battery";
    cout << replaceWords::replaceWords(dictionary, sentence) << endl;
    dictionary = {"a", "b", "c"};
    sentence = "aadsfasf absbs bbab cadsfafs";
    cout << replaceWords::replaceWords(dictionary, sentence) << endl;
}

namespace minSteps {
    int minSteps(int n) {
        // dp[i]表示打印出i个A的最少操作次数
        vector<int> dp(n + 1, INT_MAX);
        // 根据题目的描述
        dp[1] = dp[0] = 0;
        for (int i = 2; i <= n; i++) {
            for (int j = 1; j * j <= i; j++) {
                if (i % j == 0) {
                    // i = 6, j = 2
                    // 2个A加上3次复印
                    // dp[6] = min(dp[6], dp[2] + 3)
                    dp[i] = min(dp[i], dp[j] + i / j);
                    // 3个A加上2次复印
                    // dp[6] = min(dp[6], dp[3] + 2)
                    dp[i] = min(dp[i], dp[i / j] + j);
                }
            }
        }

        return dp[n];
    }
}

void minSteps_test() {
    int n;
    n = 4;
    cout << minSteps::minSteps(n) << endl;
    n = 3;
    cout << minSteps::minSteps(n) << endl;
    n = 1;
    cout << minSteps::minSteps(n) << endl;
}

namespace findDuplicateSubtrees {
    unordered_map<string, TreeNode::TreeNode *> seen;
    unordered_set<TreeNode::TreeNode *> repeat;

    string dfs(TreeNode::TreeNode *node) {
        if (!node) {
            return "";
        }
        string serial = to_string(node->val) + "(" + dfs(node->left) + ")(" + dfs(node->right) + ")";
        if (auto it = seen.find(serial); it != seen.end()) {
            repeat.insert(it->second);
        } else {
            seen[serial] = node;
        }
        return serial;
    }

    vector<TreeNode::TreeNode *> findDuplicateSubtrees(TreeNode::TreeNode *root) {
        seen.clear();
        repeat.clear();
        dfs(root);
        return {repeat.begin(), repeat.end()};
    }
}

void findDuplicateSubtrees_test() {
    vector<int> vals;
    TreeNode::TreeNode *root;
    vector<TreeNode::TreeNode *> ans;
    vals = {1, 2, 3, 4, 0, 2, 4, 0, 0, 4};
    root = create_treenode(vals);
    ans = findDuplicateSubtrees::findDuplicateSubtrees(root);
    for (auto tree : ans) {
        cout << TreeNode::print_tree(tree) << endl;
    }
    cout << "________________" << endl;

    vals = {2, 1, 1};
    root = create_treenode(vals);
    ans = findDuplicateSubtrees::findDuplicateSubtrees(root);
    for (auto tree : ans) {
        cout << TreeNode::print_tree(tree) << endl;
    }
    cout << "________________" << endl;

    vals = {2, 2, 2, 3, 0, 3, 0};
    root = create_treenode(vals);
    ans = findDuplicateSubtrees::findDuplicateSubtrees(root);
    for (auto tree : ans) {
        cout << TreeNode::print_tree(tree) << endl;
    }
    cout << "________________" << endl;
}

namespace findTarget {
    unordered_map<int, int> map;
    bool ans = false;

    void dfs(TreeNode::TreeNode *root, int k) {
        if (root == nullptr)
            return;
        if (map.find(root->val) != map.end()) {
            ans = true;
            return;
        } else {
            map[k - root->val] = root->val;
        }
        dfs(root->left, k);
        dfs(root->right, k);
    }

    bool findTarget(TreeNode::TreeNode *root, int k) {
        map.clear();
        ans = false;
        dfs(root, k);
        return ans;
    }
}

void findTarget_test() {
    vector<int> nums;
    int k;
    TreeNode::TreeNode *root;
    nums = {5, 3, 6, 2, 4, 0, 7};
    k = 9;
    root = create_treenode(nums);
    cout << findTarget::findTarget(root, k) << endl;
    nums = {5, 3, 6, 2, 4, 0, 7};
    k = 20;
    root = create_treenode(nums);
    cout << findTarget::findTarget(root, k) << endl;
}

namespace printTree {
    int calDepth(TreeNode::TreeNode *root) {
        int res = -1;
        queue<TreeNode::TreeNode *> q;
        q.push(root);
        while (!q.empty()) {
            int len = q.size();
            res++;
            while (len) {
                len--;
                auto t = q.front();
                q.pop();
                if (t->left) {
                    q.push(t->left);
                }
                if (t->right) {
                    q.push(t->right);
                }
            }
        }
        return res;
    }

    vector<vector<string>> printTree(TreeNode::TreeNode *root) {
        int height = calDepth(root);
        int m = height + 1;
        int n = (1 << (height + 1)) - 1;
        vector<vector<string>> res(m, vector<string>(n, ""));
        queue<tuple<TreeNode::TreeNode *, int, int>> q;
        q.push({root, 0, (n - 1) / 2});
        while (!q.empty()) {
            auto t = q.front();
            q.pop();
            int r = get<1>(t), c = get<2>(t);
            res[r][c] = to_string(get<0>(t)->val);
            if (get<0>(t)->left) {
                q.push({get<0>(t)->left, r + 1, c - (1 << (height - r - 1))});
            }
            if (get<0>(t)->right) {
                q.push({get<0>(t)->right, r + 1, c + (1 << (height - r - 1))});
            }
        }
        return res;
    }
}

void printTree_test() {
    vector<vector<string>> ans;
    vector<int> nums;
    TreeNode::TreeNode *root;
    nums = {1, 2};
    root = create_treenode(nums);
    ans = printTree::printTree(root);
    for (auto s : ans) {
        for (auto t : s) {
            cout << t;
        }
        cout << endl;
    }
    cout << "_______" << endl;
    nums = {1, 2, 3, 0, 4};
    root = create_treenode(nums);
    ans = printTree::printTree(root);
    for (auto s : ans) {
        for (auto t : s) {
            cout << t;
        }
        cout << endl;
    }
}

namespace isPossible {
    bool isPossible(vector<int> &nums) {
        unordered_map<int, priority_queue<int, vector<int>, greater<int>>> mp;
        for (auto &x : nums) {
            if (mp.find(x) == mp.end()) {
                mp[x] = priority_queue<int, vector<int>, greater<int>>();
            }
            if (mp.find(x - 1) != mp.end()) {
                int prevLength = mp[x - 1].top();
                mp[x - 1].pop();
                if (mp[x - 1].empty()) {
                    mp.erase(x - 1);
                }
                mp[x].push(prevLength + 1);
            } else {
                mp[x].push(1);
            }
        }
        for (auto it = mp.begin(); it != mp.end(); ++it) {
            priority_queue<int, vector<int>, greater<int>> queue = it->second;
            if (queue.top() < 3) {
                return false;
            }
        }
        return true;
    }
}

void isPossible_test() {
    vector<int> nums;
    nums = {1, 2, 3, 3, 4, 5};
    cout << isPossible::isPossible(nums) << endl;
    nums = {1, 2, 3, 3, 4, 4, 5, 5};
    cout << isPossible::isPossible(nums) << endl;
    nums = {1, 2, 3, 4, 4, 5};
    cout << isPossible::isPossible(nums) << endl;
}

namespace imageSmoother {
    vector<vector<int>> imageSmoother(vector<vector<int>> &img) {
        vector<vector<int>> ans;
        ans = img;
        int cols, rows;
        rows = img.size();
        cols = img[0].size();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int sum = 0;
                int count = 0;
                for (int k = -1; k <= 1; ++k) {
                    if (i + k < 0 || i + k >= rows) {
                        continue;
                    }
                    for (int l = -1; l <= 1; ++l) {
                        if (j + l < 0 || j + l >= cols) {
                            continue;
                        }
                        count++;
                        sum += img[i + k][j + l];
                    }
                }
                ans[i][j] = (int) (sum / count);
            }
        }
        return ans;
    }
}

void imageSmoother_test() {
    vector<vector<int>> img, ans;
    img = {{1, 1, 1},
           {1, 0, 1},
           {1, 1, 1}};
    ans = imageSmoother::imageSmoother(img);
    for (auto line : ans) {
        print_vector(line);
    }
    cout << "_________" << endl;

    img = {{100, 200, 100},
           {200, 50,  200},
           {100, 200, 100}};
    ans = imageSmoother::imageSmoother(img);
    for (auto line : ans) {
        print_vector(line);
    }
    cout << "_________" << endl;
}

namespace widthOfBinaryTree {
    int widthOfBinaryTree(TreeNode::TreeNode *root) {
        unsigned long long res = 1;
        vector<pair<TreeNode::TreeNode *, unsigned long long>> arr;
        arr.emplace_back(root, 1L);
        while (!arr.empty()) {
            vector<pair<TreeNode::TreeNode *, unsigned long long>> tmp;
            for (auto &[node, index] : arr) {
                if (node->left) {
                    tmp.emplace_back(node->left, index * 2);
                }
                if (node->right) {
                    tmp.emplace_back(node->right, index * 2 + 1);
                }
            }
            res = max(res, arr.back().second - arr[0].second + 1);
            arr = move(tmp);
        }
        return res;
    }
};

void widthOfBinaryTree_test() {
    vector<int> nums;
    TreeNode::TreeNode *root;
    nums = {1, 3, 2, 5, 3, 0, 9};
    root = create_treenode(nums);
    cout << widthOfBinaryTree::widthOfBinaryTree(root) << endl;
    nums = {1, 3, 2, 5, 0, 0, 9, 6, 0, 7};
    root = create_treenode(nums);
    cout << widthOfBinaryTree::widthOfBinaryTree(root) << endl;
    nums = {1, 3, 2, 5};
    root = create_treenode(nums);
    cout << widthOfBinaryTree::widthOfBinaryTree(root) << endl;
}

namespace strangePrinter {
    int strangePrinter(string s) {
        int n = s.length();
        vector<vector<int>> f(n, vector<int>(n));
        for (int i = n - 1; i >= 0; i--) {
            f[i][i] = 1;
            for (int j = i + 1; j < n; j++) {
                if (s[i] == s[j]) {
                    f[i][j] = f[i][j - 1];
                } else {
                    int minn = INT_MAX;
                    for (int k = i; k < j; k++) {
                        minn = min(minn, f[i][k] + f[k + 1][j]);
                    }
                    f[i][j] = minn;
                }
            }
        }
        return f[0][n - 1];
    }
}

void strangePrinter_test() {
    string s;
    s = "aaabbb";
    cout << strangePrinter::strangePrinter(s) << endl;
    s = "aba";
    cout << strangePrinter::strangePrinter(s) << endl;
}

namespace checkPossibility {
    bool checkPossibility(vector<int> &nums) {
        int n = nums.size();
        for (int i = 0; i < n - 1; ++i) {
            int x = nums[i], y = nums[i + 1];
            if (x > y) {
                nums[i] = y;
                if (is_sorted(nums.begin(), nums.end())) {
                    return true;
                }
                nums[i] = x; // 复原
                nums[i + 1] = x;
                return is_sorted(nums.begin(), nums.end());
            }
        }
        return true;
    }
}

void checkPossibility_test() {
    vector<int> nums;
    nums = {4, 2, 3};
    cout << checkPossibility::checkPossibility(nums) << endl;
    nums = {4, 2, 1};
    cout << checkPossibility::checkPossibility(nums) << endl;
    nums = {3, 4, 2, 3};
    cout << checkPossibility::checkPossibility(nums) << endl;
}

namespace constructArray {
    vector<int> constructArray(int n, int k) {
        vector<int> answer;
        for (int i = 1; i < n - k; ++i) {
            answer.push_back(i);
        }
        for (int i = n - k, j = n; i <= j; ++i, --j) {
            answer.push_back(i);
            if (i != j) {
                answer.push_back(j);
            }
        }
        return answer;
    }
}

void constructArray_test() {
    int n = 3, k = 1;
    vector<int> ans;
    ans = constructArray::constructArray(n, k);
    print_vector(ans);
    n = 3, k = 2;
    ans = constructArray::constructArray(n, k);
    print_vector(ans);
}

namespace findKthNumber668 {
    int findKthNumber(int m, int n, int k) {
        int left = 1, right = m * n;
        while (left < right) {
            int x = left + (right - left) / 2;
            int count = x / n * n;
            for (int i = x / n + 1; i <= m; ++i) {
                count += x / i;
            }
            if (count >= k) {
                right = x;
            } else {
                left = x + 1;
            }
        }
        return left;
    }
}

void findKthNumber668_test() {
    int m, n, k;
    m = 3, n = 3, k = 5;
    cout << findKthNumber668::findKthNumber(m, n, k) << endl;
    m = 2, n = 3, k = 6;
    cout << findKthNumber668::findKthNumber(m, n, k) << endl;
    m = 9895, n = 28405, k = 100787757;
    cout << findKthNumber668::findKthNumber(m, n, k) << endl;
}

namespace trimBST {
    TreeNode::TreeNode *trimBST(TreeNode::TreeNode *root, int low, int high) {
        if (root == nullptr) {
            return nullptr;
        }
        if (root->val < low) {
            return trimBST(root->right, low, high);
        } else if (root->val > high) {
            return trimBST(root->left, low, high);
        } else {
            root->left = trimBST(root->left, low, high);
            root->right = trimBST(root->right, low, high);
            return root;
        }
    }
}

void trimBST_test() {
    vector<int> nums;
    int low, high;
    TreeNode::TreeNode *root, *ans;
    nums = {1, 0, 2};
    low = 1, high = 2;
    root = create_treenode(nums, true);
    ans = trimBST::trimBST(root, low, high);
    cout << TreeNode::print_tree(ans) << endl;
    nums = {3, 0, 4, -1, 2, -1, -1, 1};
    low = 1, high = 3;
    root = create_treenode(nums, true);
    ans = trimBST::trimBST(root, low, high);
    cout << TreeNode::print_tree(ans) << endl;
}

namespace maximumSwap {
    int maximumSwap(int num) {
        string char_array = to_string(num);
        int n = char_array.size();
        int maxIdx = n - 1;
        int idx1 = -1, idx2 = -1;
        for (int i = n - 1; i >= 0; i--) {
            if (char_array[i] > char_array[maxIdx]) {
                maxIdx = i;
            } else if (char_array[i] < char_array[maxIdx]) {
                idx1 = i;
                idx2 = maxIdx;
            }
        }
        if (idx1 >= 0) {
            swap(char_array[idx1], char_array[idx2]);
            return stoi(char_array);
        } else {
            return num;
        }
    }
}

void maximumSwap_test() {
    int num;
    num = 2736;
    cout << maximumSwap::maximumSwap(num) << endl;
    num = 9973;
    cout << maximumSwap::maximumSwap(num) << endl;
}

namespace findSecondMinimumValue {
    int findSecondMinimumValue(TreeNode::TreeNode *root) {
        int ans = -1;
        int rootvalue = root->val;

        function<void(TreeNode::TreeNode *)> dfs = [&](TreeNode::TreeNode *node) {
            if (!node) {
                return;
            }
            if (ans != -1 && node->val >= ans) {
                return;
            }
            if (node->val > rootvalue) {
                ans = node->val;
            }
            dfs(node->left);
            dfs(node->right);
        };

        dfs(root);
        return ans;
    }
}

void findSecondMinimumValue_test() {
    vector<int> nums;
    TreeNode::TreeNode *root;
    nums = {2, 2, 5, 0, 0, 5, 7};
    root = create_treenode(nums);
    cout << findSecondMinimumValue::findSecondMinimumValue(root) << endl;
    nums = {2, 2, 2};
    root = create_treenode(nums);
    cout << findSecondMinimumValue::findSecondMinimumValue(root) << endl;
}

namespace flipLights {
    int flipLights(int n, int presses) {
        unordered_set<int> seen;
        for (int i = 0; i < 1 << 4; i++) {
            vector<int> pressArr(4);
            for (int j = 0; j < 4; j++) {
                pressArr[j] = (i >> j) & 1;
            }
            int sum = accumulate(pressArr.begin(), pressArr.end(), 0);
            if (sum % 2 == presses % 2 && sum <= presses) {
                int status = pressArr[0] ^pressArr[2] ^pressArr[3];
                if (n >= 2) {
                    status |= (pressArr[0] ^ pressArr[1]) << 1;
                }
                if (n >= 3) {
                    status |= (pressArr[0] ^ pressArr[2]) << 2;
                }
                if (n >= 4) {
                    status |= (pressArr[0] ^ pressArr[1] ^ pressArr[3]) << 3;
                }
                seen.emplace(status);
            }
        }
        return seen.size();
    }
}

void flipLights_test() {
    int n, presses;
    n = 1, presses = 1;
    cout << flipLights::flipLights(n, presses) << endl;
    n = 2, presses = 1;
    cout << flipLights::flipLights(n, presses) << endl;
    n = 3, presses = 1;
    cout << flipLights::flipLights(n, presses) << endl;
}

namespace findNumberOfLIS {
    int findNumberOfLIS(vector<int> &nums) {
        int n = nums.size(), maxLen = 0, ans = 0;
        vector<int> dp(n), cnt(n);
        for (int i = 0; i < n; ++i) {
            dp[i] = 1;
            cnt[i] = 1;
            for (int j = 0; j < i; ++j) {
                if (nums[i] > nums[j]) {
                    if (dp[j] + 1 > dp[i]) {
                        dp[i] = dp[j] + 1;
                        cnt[i] = cnt[j]; // 重置计数
                    } else if (dp[j] + 1 == dp[i]) {
                        cnt[i] += cnt[j];
                    }
                }
            }
            if (dp[i] > maxLen) {
                maxLen = dp[i];
                ans = cnt[i]; // 重置计数
            } else if (dp[i] == maxLen) {
                ans += cnt[i];
            }
        }
        return ans;
    }
}

void findNumberOfLIS_test() {
    vector<int> nums;
    nums = {1, 3, 5, 4, 7};
    cout << findNumberOfLIS::findNumberOfLIS(nums) << endl;
    nums = {2, 2, 2, 2, 2};
    cout << findNumberOfLIS::findNumberOfLIS(nums) << endl;
}

namespace findLengthOfLCIS {
    int findLengthOfLCIS(vector<int> &nums) {
        int ans = 0;
        int n = nums.size();
        int start = 0;
        for (int i = 0; i < n; ++i) {
            if (i > 0 && nums[i] <= nums[i - 1]) {
                start = i;
            }
            ans = max(ans, i - start + 1);
        }
        return ans;
    }
}

void findLengthOfLCIS_test() {
    vector<int> nums;
    nums = {1, 3, 5, 4, 7};
    cout << findLengthOfLCIS::findLengthOfLCIS(nums) << endl;
    nums = {2, 2, 2, 2, 2};
    cout << findLengthOfLCIS::findLengthOfLCIS(nums) << endl;
}

namespace MapSum {
    class MapSum {
    public:
        MapSum() {

        }

        void insert(string key, int val) {
            int delta = val;
            if (map.count(key)) {
                delta -= map[key];
            }
            map[key] = val;
            for (int i = 1; i <= key.size(); ++i) {
                prefixmap[key.substr(0, i)] += delta;
            }
        }

        int sum(string prefix) {
            return prefixmap[prefix];
        }

    private:
        unordered_map<string, int> map;
        unordered_map<string, int> prefixmap;
    };
}

void MapSum_test() {
    std::shared_ptr<MapSum::MapSum> map_sum = std::make_shared<MapSum::MapSum>();
    map_sum->insert("apple", 3);
    cout << map_sum->sum("ap") << endl;
    map_sum->insert("app", 2);
    cout << map_sum->sum("app") << endl;
    map_sum->insert("app", 2);
    cout << map_sum->sum("app") << endl;
    cout << "++++++++++++++++++" << endl;

}

namespace checkValidString {
    bool checkValidString(string s) {
        stack<int> leftStack;
        stack<int> asteriskStack;
        int n = s.size();

        for (int i = 0; i < n; i++) {
            char c = s[i];
            if (c == '(') {
                leftStack.push(i);
            } else if (c == '*') {
                asteriskStack.push(i);
            } else {
                if (!leftStack.empty()) {
                    leftStack.pop();
                } else if (!asteriskStack.empty()) {
                    asteriskStack.pop();
                } else {
                    return false;
                }
            }
        }

        while (!leftStack.empty() && !asteriskStack.empty()) {
            int leftIndex = leftStack.top();
            leftStack.pop();
            int asteriskIndex = asteriskStack.top();
            asteriskStack.pop();
            if (leftIndex > asteriskIndex) {
                return false;
            }
        }

        return leftStack.empty();
    }
}

void checkValidString_test() {
    string s;
    s = "(((((*(()((((*((**(((()()*)()()()*((((**)())*)*)))))))(())(()))())((*()()(((()((()*(())*(()**)()(())";
    cout << checkValidString::checkValidString(s) << endl;
    s = "((((()(()()()*()(((((*)()*(**(())))))(())()())(((())())())))))))(((((())*)))()))(()((*()*(*)))(*)()";
    cout << checkValidString::checkValidString(s) << endl;
    s = "()";
    cout << checkValidString::checkValidString(s) << endl;
    s = "(*)";
    cout << checkValidString::checkValidString(s) << endl;
    s = "(*))";
    cout << checkValidString::checkValidString(s) << endl;
}

namespace judgePoint24 {
    bool calculate(std::vector<double> numbers) {
        if (numbers.size() == 1) {
            return std::fabs(numbers[0] - 24.0) < 1e-6; // 检查是否接近24
        }

        for (size_t i = 0; i < numbers.size(); ++i) {
            for (size_t j = 0; j < numbers.size(); ++j) {
                if (i != j) {
                    std::vector<double> newNumbers;
                    for (size_t k = 0; k < numbers.size(); ++k) {
                        if (k != i && k != j) {
                            newNumbers.push_back(numbers[k]);
                        }
                    }

                    // 尝试所有操作
                    double a = numbers[i];
                    double b = numbers[j];

                    newNumbers.push_back(a + b);
                    if (calculate(newNumbers)) return true;
                    newNumbers.pop_back();

                    newNumbers.push_back(a - b);
                    if (calculate(newNumbers)) return true;
                    newNumbers.pop_back();

                    newNumbers.push_back(b - a);
                    if (calculate(newNumbers)) return true;
                    newNumbers.pop_back();

                    newNumbers.push_back(a * b);
                    if (calculate(newNumbers)) return true;
                    newNumbers.pop_back();

                    if (b != 0) {
                        newNumbers.push_back(a / b);
                        if (calculate(newNumbers)) return true;
                        newNumbers.pop_back();
                    }

                    if (a != 0) {
                        newNumbers.push_back(b / a);
                        if (calculate(newNumbers)) return true;
                        newNumbers.pop_back();
                    }
                }
            }
        }

        return false;
    }

    bool canCalculate24(std::vector<double> numbers) {
        return calculate(numbers);
    }
    static constexpr int TARGET = 24;
    static constexpr double EPSILON = 1e-6;
    static constexpr int ADD = 0, MULTIPLY = 1, SUBTRACT = 2, DIVIDE = 3;
    bool solve(vector<double> &l) {
        if (l.size() == 0) {
            return false;
        }
        if (l.size() == 1) {
            return fabs(l[0] - TARGET) < EPSILON;
        }
        int size = l.size();
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i != j) {
                    vector<double> list2 = vector<double>();
                    for (int k = 0; k < size; k++) {
                        if (k != i && k != j) {
                            list2.emplace_back(l[k]);
                        }
                    }
                    for (int k = 0; k < 4; k++) {
                        if (k < 2 && i > j) {
                            continue;
                        }
                        if (k == ADD) {
                            list2.emplace_back(l[i] + l[j]);
                        } else if (k == MULTIPLY) {
                            list2.emplace_back(l[i] * l[j]);
                        } else if (k == SUBTRACT) {
                            list2.emplace_back(l[i] - l[j]);
                        } else if (k == DIVIDE) {
                            if (fabs(l[j]) < EPSILON) {
                                continue;
                            }
                            list2.emplace_back(l[i] / l[j]);
                        }
                        if (solve(list2)) {
                            return true;
                        }
                        list2.pop_back();
                    }
                }
            }
        }
        return false;
    }
    bool judgePoint24(vector<int> &nums) {
        vector<double> l;
        for (const int &num : nums) {
            l.emplace_back(static_cast<double>(num));
        }
        return solve(l);
    }


//    bool judgePoint24(vector<int>& cards) {
//        vector<double>nums;
//        for (auto card : cards) {
//            nums.push_back((double)card);
//        }
//        return canCalculate24(nums);
//    }
}

void judgePoint24_test() {
    vector<int>cards;
    cards = {4, 1, 8, 7};
    cout << "true," << judgePoint24::judgePoint24(cards) << endl;
    cards = {1, 2, 1, 2};
    cout << "false," << judgePoint24::judgePoint24(cards) << endl;
}

namespace validPalindrome {
    bool is_delete;
    bool dfs (int left, int righ, string &s) {
        if (righ <= left) {
            return true;
        }
        if (s[left] == s[righ]) {
            left++;
            righ--;
            return dfs(left, righ, s);
        } else {
            if (!is_delete) {
                is_delete = true;
                return dfs(left+1, righ, s) || dfs(left, righ-1, s);
            } else {
                return false;
            }
        }
        return false;
    }
    bool validPalindrome(string s) {
        is_delete = false;
        int left, righ;
        left = 0, righ = s.size() - 1;
        return dfs(left, righ, s);
    }
}

void validPalindrome_test(){
    string s;
    s = "enveorysiwkzbfngqeijeynzlfuivzsbjgwrpgcawikmvsbtmxhokubhrahzpougclcfzmmwklgxfyeovygfjwdygkevohzujhztzxyfpfajlvublakbkcwfrboxrzawwmfbnxaojiwjfiwmfsjumqitxneuagjkcasiffnbidfsmfeszbjyqwlenvrepixswlsqkablzataibfoxtooerdikycftzemaoesyjjngvczuhycyydufiedzhqslekqcvcriyqpghpazmxubtekiirixiiaaumscgoxcnolmsglnxfkpzaiiwbbymaukjofbuqcluysnworqxoxabmtbnounhwfzpicimbyiuotiklhyhavdkupsgvcywzlnorvpttfoiqzrwdnbtzwbxsowcasrpifzcgtfqvxxattcgcfogpmdymzpkmeyrxodixqbvjvvbrsftcffzimetikzgmzuadutalxkuzurnnqvrkjgxxkxigmojmrzotafxhpxffblrwwqzcgqqajncmnppucaasqcwlmaxjnwmwvhlrfqbshbnaampcnrrzbiincadsicvwivbelacwqpvkkrvkukunweffgcwieiiqvuxtzbccikemybtlpnckceqdnyghfuwkkaigprxzbgqvbhwibbzabpkpnbkimowbbbfcipdnwbjuockoxshafnbaflqsulpivltcatubkgfrbryzkeixiyjrmxjkivhvozucocpcshoxpajzeftxmbtbufyzyosvrjsodktjwclrhsyywyqiywojcbazcfozhyuyqudptkdtjqfsgcdhzghnlcfroubmpclsdcuuaeazlpgoanunnxrsvbnqzfdgoasxljsnxmpqhsegatrvfwgncdzdgvibicjbhnprihrvxfakdegoexfhvnxqhhpnjmcvbatrllupsyijamavydakqknftfsjorrtgkykbkmgcrerwwgzolbxdtcpthkxspkamvkaklzpticptgiicakcmylcrbuuxzenmcozmlqjsnnojqssclrcofkjxwxcmlvieibvigmnffhfmdyedtdbfeuggnjcupdtlfnxfevxophdjhdvzcwsbernmoorlgnnohlapsfzkntksxpcvwotetbpbixmibnbidrzegygzqzchwvqhrhctzktsxnuyfvfdcsxhjsbczaipggesjnpedbswgsgnmcmhuyjcimdfqlxvwjakiuninralqdmorvifzlgziwgcocputjujhclyxzzjnbhxdvkezjnjykhdvgbnvitegnqnqyqxtxivuupwhqfncgopiddgsexdiirlcffkkqtytrkbzpetvbysknlilvnzkguvocsgwailskyowoglnnpjhdporzklwasiaivbyzutnqnbfhtamxqawkqwdfirygvmszzxfjlmbnlaiatjkmkzizrysnzeffefwejpntamjgboaplakxlfmupihwxmvsjlhokccfdvfikdmdhacakhhsaaspmdsayoixrznxzqosyklnucghejhxubbiznyriuxnfocrnbcfkvewxmdusqncplxtmnpqjohxybcppucekrrptviudavgxvilekhvsqajstazzqwevtlrnzwcmgdqzqxzaogpshrisuucrkxthujfzscmmszdvoyvnawnqtddabbmxzepjyzrnyslgncptblotchkzjtumpmbwevnfnolcnuryfaugbxwsaqlevsyptxdhgbsyoizyuavjapaobxpezujumyiwisvmocjlfjlfimviszjsdskigceygrvrmcbthrghlnlbduhtrgwedzjqcpuxlugatlfsvqzppmnnrqxrubwwrjpfiugenwkqhjkcsfeqgjenbhntajsdglqfgflprcryfquqkutyltkwuejoeziwgpzjikfiywytihgqzxvqwrokvbcwosuipxkwxawnptbndvtvdrphlfwkxfthlpamhzaexescnttqpuwdzuitfiyfzzvxqklplfrimyaoulzexzcjnzgwrcmqwzlyvwuxssufjjlhwtjrtxvkqyvaktmhzkrsjiubdjsjazajthxfiwvfytjnpauuzjhkdkxcyqgpzhlpyrmlykihlrgjdyxtvrraoamivkjwgjqamxbhbstbifueiuzmjjncuysrmfivllrygrolwfocfdpzhwufvirdifnubznrolrbdoofmhubqtrkmydyrjnaykcdsguhrdjlseoashlriojinuayfurqqwzvhgohfrqdfbrgszpwqycscebqigrvurrcqsbyvkrrmtjxdblsszdnbgwjzjzvrcwxkdksboszcnxpoznrrkxsgkerhihckqcqkzunfkzmlhyzsgmnorlblorngxnrlgdilsrhvzxcdyzyjzixtkuivocftguzifpijjuuyjsxetxyjzhczbarnjzwhtkhtcwznltwygndvtcrieuimeemskgbjhygizgmpkawkmjgbywsapzufwcmrmlvjtbqpugchvefajwholukgahhmvelqmvszcfqupfvrergtuqcnwcfyeatlnpoaknvagwswsljpbtqgehnwjtrrwupshsgcmeaiwoqfrqimxhaynqunapczzlmoqoqkfbedomoqvxnxpsnmerfocrildafqxldsjedtzjlovndubmgorlfitorilalklxhkonpcakpyurfdnuordyximlmeayvopxqqxcqbseztfyulwcjpboajzyrtmtaricgxezfyuydjjfvxjkzfdodbgezcxwmcguaxslbetbdlzvymmodssxacwtsdokniqlkcfeyhbewwiyatwiuptrdggrzcfxxuoapgqpjfpjdrdjnvnswywfcqivbmpxfmmdriyevzgfrhpzpyajjwukkkxqkfvmwtvygimdbgojzzljjqdckytntlhbwscjbjovdrmqpdnrgmvksgskvahzmvcrfhixxqwhuscihmvhtkzkqxnqttiigfyjqjppoajmlaqrhsusywkwffnxnehmatywjhnbqybwswlshyphslqvawutebuvkbcleitwdvrdkrvagjlmmhppnnhsouzsnlkwbwcvxehswwhupstaahctyhmplfafosjsxmscmkaihyvxsnwnrsnqxyicufbfbeklfaoqiaqmzommvcafxazazbzcdjiamtqomquymhtlyjhvnfwmxreuvwobqznusiylirpkxqudbbxsctbhtlrcizzqkkxsavwkqxqgociyjqwhcmwiewctwcfakhzyzpqcmtjgmckpgdstawlfzvuwniuiyxaobfkzreplcutgxtvehjruvxstsozqyopznpugdrymapyyntyabqhusrvmmdpwvxbtcghowkdsjjjywtvabhrloncttptlmuqyhzyfymcbpfzraladlfislgkxwevzqwcozmqreujpnltvxvsrbcgidslrjuqvsuhlforttbwetfexcrldreedtsroarcsejtkdpszwwzvfgthlfgcmkslcghdkbcypjudtcblowjsdhuaiskmmixbakukqebwbuxvaeulqgvkxcerlronrsnchjpvqthkmtgruqinzaipfbuhfxfjtkvzzoaaulcpikqbbydobnqdpnavfgbdbiagidjkyuepknfgyctridqncktzafxcgilbbeoshqxoasfuagfdzabajizatmjhgbyqcezrttjowkpuuppaezvdtzeccfnhbqmrnxkqlnhginrdfcakfctnxuhjzljqyysptfpgbedaxfjdkzucykvgieslcbfhrbicrtrqmvldyociglbkjkqqtzxjojxaaliusbfskobtshrgkfdjbomkishvrkmvhkyxjewezfabmzgzpofsolthducegrouketgrxquvtikovxwkowpeemjawkurjpvokndbvxyhifsaictmhwxtyigsohekhlfhcuxbzezbaeetecphepnyhhimlvdpwkvvfyttkklojbdqbedfsvxwjympdcvpwuzolsztosxklgtqokojpauqfaoojgfnstuweyvgnnmpgsmstragbxhrpwddtfpfnwkqqipcoataxueyawnxgdzuswtzcsfizrtrskjhoeeemipnbtwpntuedqrqvrprmpljsddchoygoicnwnedpqscbwecmfijtltvuthemmgabquhyzhnfxcacnzhgyzstigczpsplovxejnnhfdagxzjbhqcjlrsnuihjcomvifyqhgyscpnrcbmuossssbmggfkprbcqpkkcfsddhwgwhzpyyhqsjakzkeeykslbrkbbjcjmcdgxsmkshqisavktwwqgfzummkyjmqkqhkdjkpohyqmftveecrwsazvvaldoilnsidcahuvzghclfviwsmentlgcsbthhlsdviylpqzhcgnbqjkjixpidaogwzbaceobvdqttlhhusvlpkmovgkuymlwcibqywmjbbpztppkgwirtpsdawuabhxbxpmxggaeltxfyeszckldqtmsysjkfameokbbtwezspepaynwelnzwfiurigwjmdiumylrjxvhbrwzfdgrigljjssibboazwgazhhforookkckdihgqglacehweovnwzwyolbldimyyrocgwyyuomxdxjckpdxxpwirmvkyrhjbmvkemvnwydyzhnoqpexnijagvdkgkbrrhsxjtcuniffjruigfwajaxgbtzxlvhamvjlazpcnwfasrbpmipimomxgcnyputacmzpnkkdepfcajpcwwkijkulcvtorkltitigxrnqqhudvnxhqnheqivcicotszuatvhcalquqikirplbfckfabvygxifddxisziqmkpwhkvefmxkxqpckryumafbuqgwnceeapdakurdtctgtyyrveoormdpocxugyzlinrtfrxveukxqyybaifcshjgaaxujsphfwpgfjciehnzsbhvmpurjlzfiycnbmlqvhtvvnjsqtvcpqnkbqrcddqyoqprrvutupzhcvscvvgpcofircifsfceutpscwcfmkyfftwxgfxayfiswdmacexafhunforfzaxftmeagqapcahbarpgwmmhzufeobsnlglbewptexejpcxsupliymkllhclgqsflclrbxqgmfsvzfaravapkvgmbtazouqmbqphnfmmdqebzvzlkbfscsjprjthmmvbgvwyfhnefmjdlyzduzzcoydgxcawyxwckymjaufufeukmghnjzxejlxcqrqxzmrvhhfvepimgzxtirxrarsvimlohuzvwricbwivvcnhntrqrqubdlszcjfmdyulfnfezjgicmjvwamqvlhcqkdnafczwdmjbfiwbexbnksnfwgmfltfpocdmypfjnfuppkhahjcomjyjhqwlykggbefonvimkyyzdqthpzfrzzuzztdxktnpdyaxosdnwothxybsdeseomfbalctslqfefdrqbaiqbptqtvwrgsfzapernxlnsfuoytbzeoaxnikeblnhwibkaujtkdqxqssfslqzkbvucrfgbyksvnkonfoewzxumtmpwfqitlasraxrjqxzmqvxxknguxkwxcsmtmiocibkfuxiibjdgccykwqjygmrqbvslhdsadrhecyaddzmcrskesdgsmthqytxjbhmoeibhaoipscokfbkrunyxbgabrhtazpglphttuoqxsnteahchddiavrqglyaqezxebgnsrcxdmugnxeorobpfhuadhvuuuaudznmrfrblixirkabypmcormdvgnoemhoyeztyyzjkeohpuswjuzdsoqwmxgpkhqoibppiobjyixfjwryivslohpbjuilxdcabgymkxfxezfnvetqdhpahpfjbepuiaekmrxsiepocinskqtecdvwsnlmuxwggrnkehddgqumclikvyypgbivvxavqojnjutujsytnboxlithbcalrhqgkxhknibbignxtrpwcjeibjnvtzfbszxsbkvucdkszqwsmjuqqcvretdzdvpzajnzqkfwvmrcibelnnxlqxoebgrytbavyroiqirvzhqpsuslpqrazhxaunijjlkmiaitcqanearsyieqkuvzacrfcmvfnscshpclcwzvlmqdjwkzxagiruilgtvmvkinjlvlgoatrqbdqrtrmggvrkycojfsrtrwrscdixisrxzhwqjmdwdwyrztuzhdcspcbscmmwnzfzcxeblfeipoluhidttjucsvxhdkjrelglyhynhrxfcdamfskhyxrjoqewlizrcqulhamedbjgmfwgasqyyebdyhmdpydvduvymtidjtnurizhmivfblafyrydwasjwayuoorfwowsxwvhkziaxjolsaiaoqjrkqyvciurglhpqwfbdghzseapiueyspnzmbhxwicpxzkuosmgsaopoplsbjnbgudcwuzqlrlwikpecdzrefifxruuopyonjcwnhndodcjavujxusfbyisuuxmzouzbyapqyvgizgtbdewqjnkztorqodsqruabxgtmaptffrmommqtrmphdeztkckpuwkxcespddexpplfspxtszhklndgwzjfmljifmtymzlkmosksajkrierewydmnaojhvksqxnhctsjcwecpamzkotcdqqnahqqhxpcpsqzuntktomsbkivdmyshvachhylwkwwhbnzhqikymujfqqcdjxolmwzpoqnkrnpvdkgarvlocqazjietgltwbnitbguxicqyedgiiyfycinvqxrfsvdtfohubffjvalcnfbuusivmmeahmntaoknztejwnhubyuifxsspcplbyqlgceaivehhsgxytchnmsonyyjgebiaokhalcdmwjgqjuczrrnulpdpphvjmcwwsbvascsxpsbueetwmznfkmrgtkkugvlvaoedqhekvdfqkpxsnggopuwcdgocvbajickfwjvbffakccxwfxcnkahjirqsyguabwumcuvdgwtamlrqguraranzmwetrvhiyonnkvilpoaitkkrcxynsnydjeocphxzlvszqczbjwwpmdfdqhrirvxfdtminryvueseluzrdtciqigjglqvjlnqcbcpoomqksyxdptokelyyxjymazpfmspuaqhbzmyxxzjlydqeckmfnalsjcgbrrlaymvzpdrmwcdhkzkkqlasbroikbmqqphujytcdbryakwkkkggnxcaptffgyjxfqbhdkysmqiuzamcnsyqmpqfmnteudevkegmqygdkqmsbbqovdvqquuzarmbeisdicmepbvfwuezmogjgwksucrusksgmeykdplckfpuiyysaggwuruzqzzzjlgbhrztkihgxybflhwhwitrgqvueysvamtnugkeyovsazimimwzfnajtdwypccwufdqjacgpxdkdonndsgmdlgceegwrdubspiupyqttrjonfhgaoogylhldenopzqgfoxqzefsvtunevesraaqqckeykmxkkjmsjlyofowwmwxhivyzzyjbofledjdyhrcppdckbmnmrkbyjegheoavdfldmoigzkgqrzsnxxtqpdsrqacwhpwwqwdjehpangwinaohohyrnhdwjrkjvpetompgzfxgocjbyqgwzapbiulafimqfzxtauiwgpbxwilktdqvcrpnurfextzzouhglhhkhssnsxethxmqmvqpfmnaqlgtxuyvfdjrkzallgtaepzsowqbaronqysupkvnjnpxwclsgmlppftmvhbkmogjwekjlfhxfdodpdpkavantecrshxbtiohzccqvwgnxsdkooswlnvfdbnxeuxpizlaztgyseuguyhxbhfgjeqkhwfndtcyukmnoxpahjupbooqelwwjdndidqnnuhdmsuglwtmotaznsykdprsgezcxyqqhugfvoipbmzdfwwbbvqcjfkdnygzdhypuyavzpegoruddxvnbfwvybpdcwqlkctgxcggvgabtlsqzccpcylhhfzgajuozkhihvplwsvmwxrtgdcnognvnflweyuoxaxbpoyudeyjyapdeoralrteyzcfimglouvzvolpmdmaeqsthrtwoenaqmpvejthzqkhvhgljfjcrmchsykxlsctrizujpruikjbkwcaaoulntmkpkuenbdqlgitpyipirsduhrcjoorjvkhrmmmpwehgogwdxcevnyyvnxzepitrageixptwphdprrwxpjljthermkjhqejpdostmpsteeacetvqqcpnqitojqqnttwiinncqzfqrwqwjtddbguecgpxuxrkinindhtlzhykyayzsgjabktqmibgsfaindobvhzaodxxphgpoardpcotrsjkmzeqmzoocztmzauclegxizjgpwprbfknvxfyqzahnkwdbyvenuloplcosowxqgqvmffzkciqeckmaxzixsyskulymzgunehbatibhgddurmtzuizgqvmkfuyksefyiudfzaxzydzjxwrngdbltprqdxxvxsqcyojkokcslqnucvgxgiptxqgvpzypblfcpzouarppgkcsqjeuslnnasevlfwbgwxcdrhpveylpekbyqnwuiahpahafnibfqkgppoghlejwlxbnloentimdhuptleycbdwfcttxkipgcygcpieayfwmtibylxiaavwkxbufguxhqfculqwdewzfrxafnekxqmqbpgtdzitqdhcwaacbirxtcirniphjpwgqroxfimoziujwtscagoyszujvsonzfvzmwgfgwlmrsuauensoiukyeviqbdcfvjdhrxyghxnskxzmvlcbutesvqthlqkuamswhbkbszzkhocygijaqrrecjqdyyqdqnejfogexvahncrppumuluwamhfofpqwqodvbtgmlowbrnhezmxitcupsdwxfttcmrftrddmbwzohluxuvcjiznivwpiypddyravzaiigwmuamfverczhqvkwfotnrsyupooioufpnclmjmducyzoqeldirpzmxeevbmxzrdqewvxxexisvpvkyvxqaerrkduzladnfhcapwrgrnvvqkuswmdnwzibkkljlwevcldipkimosawgwqfrlszaciuazdryrfunfmqunxvoqkuormujcxkwxtdrdtuivwpivqgrxzifdqbfdffcdcuarfcpaquocsorawuhaijmvtgkffnmiepqzkzwkdyeymujnxpbjsmbiplghjvllunaoapkcklfocnygiuloaymmwbmakzsscrugtiomnzrjakocchlijbtkbjtokqgudxvqlydbxxarkbtrmqbgovjknfvxbxjtoaukijdongljwtgzgtdlpgwaccbacqkhzmepfwuruqkxngqgxxbzzzyfgvvhdcpinttrrbgnytmfcrmxkamvwjzyoqsetlrfsgbqnalebkdyturseqjhfyvgxnfwrtdaaacxmcbeynhvwnolxqcjxkcgjgncpalqibjkyedmkmuyyxlbncgxjvmzgrniqdcpzprdlcsiedfaqihheqtljdvfbsqcyeglasnhsabjurskfpyxakeijiwwkctjywtawanixwxtjnrlbafyglbustmbwybstgnphmtflugjrsoxtrgatgtgemabbtvxetqutgcopmsksdisthagrbyzgvdsvurkdrspnfxeptqmwfvfmprkpsvqpkqgmcnsigmncdpogmvfxgxzxeyulgnilexaootamdrcnvqirhqpcomnnzkjeztpmxqwhsidkwsyoktzapvloxzjhludftsglwtqqlruvgkkskkkqybyywqslncvlmqlmhohxeistseewveahnoemoinkrjjwaqxxznikrbdunpxncjxtrpcxzovsvrthdgnfhuhicylmehcctjzqmhkzwmphhczejkvwidacdcmvpdgoyejybdxludqpfujqtizhgalrwmghmwspqdszvzvgnxneikkozpkniungshkzkzsrvhqsamhuxkkalahonbkcgaoktjlzdmmvljjjrfipmglxniviyviwivnnbkrbpfshazhhitzksgfdbckppkkcholtxmslcixnhkdnzgonnvborudlywxbknzikegxblamwntyswkftitqpsqtjpewvcuwacybqgqhqchnbvnqgwqcyygvhrrsutisrvgfkwxoiptdlgqhpfyjaqvyybkfpsdzlnuncjeearxepqxqaptaqgdjzpqegywzblqlzsbkpeswkurucuyzalkohnkozkoohingxxpcrdqflbuqddtucqkmzslodfzzjwpfglhkfzbplydphyhqbwdchrffvhopnmwdvdlotncgjjfrdggdkblqvqiohuontyjrgcarvtkinpnfwrttlpigsdqusuaejhruhdbulrgonribknxxmqzkxttvbcfpupjjxdutiuhqyvjhpgqupwotjanzhhrbtbciskvsqzuftknsydrnokbuzehnimnrraclselgsfimkktmevhxfayazkassexuqkxrltcngglrmjoktfdaougpajwpcqzdzufhcrvxiwuwtlhyyfwncvrknhacjixwnsrzxklprspquxqclbmnnspwzghfjrabshfldxqihjxcqiltiszmurjwbyptfmfwvtrrnysadzkkzdltxudmkxsyovbchwnjbglazhrzypllbgfrkjkauftbliyndqejkrimztkcqfavjysmthbmsioiehoreauekwmcgsldzjlocrfouhwmfaghjudihgqpobvwpzfilsnprrxzjjhrbikydvsrrchhevsoxowbctcaetqaoeskefzxhvvnsedblqmyqzqshuxgjmanzcilbmsrhwkhhczxpnhblqcdwyahkodwaqamwtlelchedgxunmwcuonevlduvqlwkkoxltodkxkmkqyatwxadqxzwgluymhqqoglsyudipvtqcfsetlqxmelteegegrgdxqucmdelnotpqakypzkxwxwwzmuiamolndpftwwrtjqiulxfjnelaszkgbjwovgmknrqijugyoaaezvdlpzplhnsnssuozkuomrqvabdppksriaesgkwdypfnsfsdpkokyiqnygdfnehrfnixvjiucgndnsclursgdethzwceupdclbhxastioxzcykxdesmicfbqmcsscqdqprsawjaqluozjdrpbzbyhgttnqxzgccxgtusrjidutrrndtxxaopxrpmhxsjudqlexngyrsunqjmesvoxlknhokhtcyhjylhryaxqyudauwgrlcmmchnggjstivbqnqcrhdjxvtctjwoqwanzjqhkxbvhuckceogzienxpwyfkksuaarndiqpzsesbectesifhrmiudvxhaujojbbrdzekdyytvxguivrtntndyuwhuuwsaezekhjzxmtcetbdhzijdmbpcginbnpuxdouxohwyqoydfdkkvpxehmcffobwdbsmnxgowbjoespgsildzgdxgueeipbvbrvbzkygbabolipwiomjbrrzzphoqjihddasdjvwdpkohjkgcsgqpzlhayvduffpmyvkcyxfqmsetdzbgkeqlqgrfsnpyfdxooinwoigwymgpyqmyolerrfiedopfxrbiavedzkoayzplzjboqfmkfdzhnawjhqklrdrfcstredfuulnaujdgixcrlgnscsezrbowcuqvmqgjqepbybeafqjfhnlvcethxshyrinmtucyhdiywygjfybkdigbgkprninimgsqehsmgnwlcvmiqgfhebqmygmyquuzmgphsfxrfjytjnymixsoqvtcmccdklqsmwhxfswmphfipgpkfzgetvwtkjkbbufpoqquwjvroetvnanmqkmqvdfhnydfyairdnhxvmicsjkejtczsgjxhivvpwfnakzdlivgjyytpnjnuzqgkxmvovjgbotnlexihjtfquyuqsrqjatqdhtohfspygvezblypakdolmrzdxgrkrltwposhuerialtimmjpamhawpcseucqykfyofcygsuxrlqbjgeaqcxvybuxkduuqlkzggujpohdzcipvomfjpyezrvspyjirnuhuowjbdazgxmsxeebxrkdcfhwuxcsuryzwlayowyyivmzvjhyyrzrpjsbkyguwgvzfhbldxoaupgkinuaqfacbdkotpwseflhferxqwsaxwagyrsamuzlbztnukaarwlrbmaicpcgnzipoiaowgkbufzfgimcbormgeuvdwhccwsfwfxhbhfjurgbehtutxfbghocswswitshevmzpthgkmmhyckcrvqnijzfkdaaxqjeydmazirwgblctntuzumwzslcxjlpwkmwdrtyueauqvoqkzzaeccxrdaharvozjhtmddttljhwylnbdhlltjsuzegpdjipzyrzruynunijrggpkhgcjlpoukueshrcveufpbmtoqszucocizdgxqsmhmppumhhloajroqqieocqhownzazlcncotokfgemtxksikoonaiqdjqpqbupwbkbhogymcqwoodmgnbasggcxjleztihrgjeqrniwnopdywrzmcptddmbajslnfjtnwccitdnkbfejjycwcdufxtdiddhuohhwuncriwlabuygrrgmdxhntpahztiahbufzvjqtxxbdcnmbbpevwqsocfjvlmdveoqlgfnqxlsjgbqgbzgxezcvfufseibatauyrgafbauqpmsgtuhnuzocbpedufmebccrvnfdwxwfojbbkakfohzfvnnqyqgjludgvefsvttvytlbozpicdwrtjlpzvhhgmfzlwrhicuiybhyezrnipffixkfocsnmuwhpaswrirbvxcyqpblpxancbrwdqivremjpaunbrreepeizbednfilgajkunaovvlewqzdpcxeygmvbvurwfzletruiqmwvqmuasynenqqufmleiasyqthyrsqaiiqkubwbslhjnuoxujzrgsvqrdujxgsladwfltsduxrjwrdhrkallcignmtvkiuqevhbwcwayiqnunhekllpzdrgqnrmlnlimswlwmjwaekopndggsropatvuwdhgsmlfbeszafgqydibfdxjefelxyfbpyxtswzbwgczjedgrrqwyceciuehyomlsggksjlmwlszqhvymmsnwavwxmpxeckvsgzksejrumrkgvrtrprydzkpqpxpaxaeepdaonkeshgvjqqiizayajlpnruwbtofsiokckqqhlizccllhmvfvupzheoiykfkmjoyptjouiwmvlxguvjjhzcdosmaoevgfavkovzynctihkkholqrjgeeojhpguyrrshdaydjzagyzipjrylsbxblungtnvjbuxlkvjzappgtuvgwnbryuspcfziwnfvmcurbjgaeoobmjqokzbvmscxfzakamokdpqpergpzsyxdyhfhwdgbewpxkbbnggoomawgakhjxohhxtwvoiaujlzjzauuicshyjkdkkeuhirimvtthjfbsqahpkuvztbqyivyefqlgfajuuxqkttusnfcjynaavaoyceevhfrbmrythsfkltrxtbnquntphkyukmociepslcpirzfclnqjwrzwrquglqhscroxhgrcuxrsaskvcmhvzsdgavvzhmlylnjscpulkdxpbtrhtllrrnpinyhpqtycjavoiurnqdwxyrjkphtmkwycxzmdwmclilzurzvupxiasczoxexxwwdfpmtelswnaajknzehoiqvismvszovvkxhotdpbhiborroqqymqblebdhjdxsvqosvzecvbovcbfdcefvtczsbagyaqwweplsadhmtnqznsyotcuhyvkwrzcaxjtirgaqihwzsvbahfbvruvkbvvqpbxfbmwjjpxdomcjriikzwdpgwsmpnhopetaubqrwfjedxnqhruknotjrcufrlzviwcxiykidsuwutufmhhxxojqmzibaqjrnhcerjzapborxaqeabugsxtngqugamjsmfqepwkciwzioprptmsrcwcjtmzlfmuvhktfbmzmdfbcizpoqongajwuxjkljkibubewuroufpcqwmfkoythivhyicjcyutsmfegmnotkmdpzsfqxaglweplcmhsxihkcxavtkionryfvrmvsfvtonjzlsmzvlfedsgqxuqmzgmwakmudkmbdpfritfdntycbzeeqtjauzsmjzbcpybrqgizsjjmanpctequdimbadewceqoynrempvekfzedazewjqrluijsageycwmzhhzqmcusediqjyqkclpuxgtykijjntkysbdbegeygoxgubbxdmjtmhkwjosxtusosgqjdxdafpbenlwgbwiovnhvkohumszxbdrxrqncxwxwrcozanlbzfnkdydziuqxlgonmchszmmpddghcbxjdxtebbnmuemhirdwzodyzrougmxcdzcogbkgyskidksegrxhvnapkebfacnaaxeqogffvjkspvvbdyhokdtvfwsiquaoukfykeddiggfkgqnkwlwzgflyrqxkfkkhzbdpkfajetythlsjpjizbsvgozazsvxiqgucibrxyjdhygywgcmobqptinbgwlkjlmkqwypbonqmcxzkuedylmedxolpqbjbimniajattdxtzcgzudrxlbdbnfsxvwcxxhknbzjvgotspiirdpeanuwyfmnipexspqhmgjgojwzmimqndnaqkwsqadlhwkndfodhzkfjduwxggnpbewgwegwwrzissatsflmnyxitpzpcnszkrewazuzmpuvldreataesuhvdtvpnekjckgjzdqgqgtarhqvxmioogfsmfxanuaavvzevxnduxzulagobgrcpiryrnqzgvcvqyzoypzhdwwyznmdiglhjmmuloziuvtmmeikvyokwlwncjehzzpwnxcecrfxmwmqvzinhggpvexlrpjqerlbattpdvddvinsjarizmjpinsktbdxtbcxjdjztksuwhuyqechuwteiswhtyifxrsetmryfkqvotsgcdclvoqjikfbjgtwbdjtkhddboieyuvqyqwhjuqdciejljonteotkbfobcncryfqfyymlkdcmymtsfzpuzfxccmtygzurslrfrutwspwamblnzltbnsfzbekuxwlaebbvsqpkhsqndfqjgvsefnavvehcbiwmebgidmihjiiuabnpihzrkiqlumcrguijvodewkrcnpicefaqqvniciowtolrxhakdxzggtpyygibomlldjvaogbnxdfhpniisjtjkuceguzrhhzgmffswoajriqvyhcbuuxfqgbbwdejyzhqgunanlhoehggtnqencuzbmweyfahfffzwyelalamyzvmhfnreitwilggyiycgjkphgyuwctrpajqciaosarxgytiawzhhxczaoxwtxkelbnpbmjrwkiyubdvjcxudmxtyheqlcgtpwtpupotzvktupqbxjvvxqwmcrgzafaekvskfgmkayqyjmcelkqroluqjcpckstundjuajncvbjrtyrwmrkdunwijqbsersxxgcqpcaetvrhbgtreplhixmgnghatsfwscrnnuaiuesuokyefckgyxpdiyrzedmaiotpojhkgfrbjzlvnpiftyrjmwktrbtflogngngdvkzwqtkwltabyefhgsouztrxtacfxchamilwooajcvfyjokxsdyuicnsohxdscaieavptvigbkedwnvmdnzyoohewfkpwzjzooquwmqrgvshkxvvlamqqiongyxqvrqfzasciwwzexptfhsmhfqzpkoeanikbzhcaaabbmgptdkbaxdxuhgnsuvqbsmexftoavzslfgbvlpwgxrptgmnjvbbbqfdbmsnwvrudtpykhbrxgrbglucdpycnjitqoeelbhixhdnhttetksqavkjotfhfmqdnccckzidvimopzztrjkldyqhzdaiigrbltynbwauvtmhbapykbqfrgqykriumotqxmzlcbsmiyrqvyfzdcoutuycnyuskbqvznnirauxbxgjgbtowpiuywupdxbqbecdcfjtzfmttprenovbylhugztxvqbdhvkbiltqptcnblibpsbatrisivrzkmutwbpivssiqgxzdpjctfsgcabfgsdsotbxctbstzgpytjbyclpzgmyoiwujfumtmpgmxuhowtrwpymuvybjcwsodnqogpkdkxxeitocgonphdvkomzamydnibkrancupixjmlydnykwrehinslaxmnwqoukltvqfydvzntmcheeqavdztoksxwkhybumtgzrnzfbmsbosvafikqtktfpnvaatdudgevczyyxgyspkawlybmimfohbevkvsnphmdqfrxnfvbrcxjrbyxvynjyvoavwvgrsdbzzwybdavewylrandfejstqbldhrzzhsegsxxjxmsxrypoppjycyzmpbogiirjdlldlexzxnhlhgehuctspmbevykcjgvfheqelyanfrbkzngqvswqcrgshjgogwqudphilstjtvodoivkfjauiduxuhpotoohgyotsyukoensfdoqqeerxuguwbyrugpewgfjgwhcrfeblcognnxofjkxvmtqwrkihgbwdmhgvjlfigpjchvrwjcmmorblgqhvrqnxgghaufhkoomxqzhbwhmxswxzdgqdmjjoaqssjltsvxpzujpusivtuvsnqogeeiyshfjahyvgfsntffkghhwixphcbfvhvoqlphpupcjyqecsnchsejbcbzbibiwsmdfvwhdfbidosrgltzuqtkhibsjghuvzqhkrkjszllmehtoumlrniciukrlgqkgvbosqwjwqnztmzxmqewyxrprcmawvjenxhystywxvymtczpfqteopphvxpqwifmipyvsjrrlnavzffkempdhhdtktqsskbadpruvvjggueoykzthiazonwdfscciemsxicwbaireemcyxbfcqhotoqkrvfosgjhvgzribvjrnepogsxushqeggminfoibkrmcivyqzgushlkepudszbnffctllxuaqjfzqpngvoqdoftnfglzhajavbwyudielgakspsnwlbsomvcgybqcyyffbqeoqsoursksfeapohitedgousowsmtlsxutxpanzxdcwwyjpirwrmbtwzgfwjmwdffazuomvodzlhmevpkpdpsjqyfpyjhqvxhxfhoftafjayuyqwugiyrqzvgiwzpwapqwjzkwvkojsrlsbmxvaixfkziqzbftmtkcvzagvenagiuptpvrbjyabemnjlsfcyudduhpoogvrgohzkxekmtnnsthhaqxgwegexijktnnvtivweyxfvysvxtrwjrmdtgtuxwjtripxepfbkczorrviqklhunwfspvgtbrmbxescphwngznkeoqlkqggflbxazjbcbfgnjetmftrzsnagbxqbquqzwogbadrqqvdkdeeqvolrjuhmaloevblkoktdpfzrncnwhscijyxyadjtumofhcxcjcxqdxoclqjmayktodfxglrvaxkleuenkyghshxlizyxgneowteddcpjojaeenfgxalrvynvfkneeojtvwydeuntyktjtnwytvlzfoubbkvoondwcbskqsiuuxkajeufmfvyezjrclprnlfmmbewpwckozojwqzokbypggahutcjfoosqowdgpfocpupylneusrhatepkqmvcsoeqsmqndoonowccccwwqcxtbvtygermbpjhcrjmsajrnfmmobcdslzfhaxiaufpsnijbajrrfkewmqhsklsjtnmriuqlhmxwbfnpzisthpdrumtuuhhbyvedbongxtjntncynslwxescygququogpsbgognhoudxqsqnfstinehmddbqhkfgcwyrbltmhnxedoueykohryzeqkdwwppmnecbrqyjdmwrwtlgtgigxhtkyxqgbqanaecsfgwyjxkoqotbtuogonxfxuhhgvmmplpixwjpsikchdtetcryzhfwodmzctzfdadljbsoxlaaoynslvvegzkdmczvmfrwbcxtgktoiswklpcxmljflvlpyjrbfufwsfjfhtpntgkeoohnnuoctlvikukxvbtybeecqgjktoxipneostaaoquoxzldwycoylpdbxgxmfgsohbykaoyoakulvqvvsvxwscnuubgfczoavjuyafxqeydzslsgnpqxkgnwtldlhwfgwqvnfakwljsorooivozhsxcezgdaihyrxpisxydyjfsckvceflholyddzycogmjykhqvsgohvelgpklqgciyeqyijllrmiwumksdmjftwadsbeeqjydjggxmcospaeoshoqivwtcveoxundpqbbzbrtqypkakolpbnpbpwcykmugzbppijqnnianqwytgpstbrrfikrnkjicmqbfeqnddmsvdsarkftqxbvscndjvixvhlfsacxxxoxkpvodpysukpqwubnntdjlqraadvkkqpxffeqhcmsyxrevscyqbbjuwpqrqeaahutscdicklxuqqmmxsrfwlkipxhuamdrwnaerrpdoceavegqyqvhevqvuobjoyyvpnphnjlizokvberelgmztwhinlsavjmsukwlkwysrzzimwkltnppntkbqqynlcstdtlrgerqcivjsappzaycowxhbifdxsusubfvbrrkixvskdvieyoearrccscrfnktztugeridvzhuwijufculcwbvlislksnwiamcbcgnroercnaraebwzdoeyweqarzsskuxipixbyexlxqldskiyjjhfauratcmvlrbqbiforueoxckqhueodijoccjeupyjinofszspurjuazgadhbiucmmetxysodhabepvnmgwywyusuccaeneydpwhkbvyxvcqfdnklgyjamqgfyidomnvbobokbwvdqhccacmeqkdquyagefeqxoylnnbgigvqtddqlfuishhssyrsfkhxltiqkxfnpocywtqidlhexvuwoyjdstbvysbntvbijncwhcaruyfxpuwagewlxouvaddrzboumpbuxjtmshrfkfprhtwydjeeoezbhgckbqvyegjoqthspihwfkyoxlthnthxkdjzhvupxyqziqzackhqtnfvemiejevwnxcmhtjjipiyhgzfoevxkduisooptkvloxgbjhdcqxleucnkzbzsrjayusjgiafegcqhlljwsudcrmmapqhulzivstyexzdlfpiaplvzqpunvqewahrhiptlnjbcizvpjuqhsumgebahaokkivfuykuywtjejeyzyhbwtgwlwohbqaaqmovmjmevxncdkqssoetaouizhlmxkrpnghmigswbkmjxmfjmqtwxwmoepvffxzjopcrscqzgjpwjkbdddhcmkcmjzznyokfcaticlltvoczoqqvaejteymjakvyrkqjdovksqxzzqipvqqauuuelthnkntuwvngdnslngznfnajevpyiqsmbbxwsihbqkkcquokygtlehvuaebaxfwjujfepjhjriptbundwldijanuygvixjtxjbqbalhjoszcxcyqssdvqgmcxkoejcdhcvuaujnzhvrwueyenlsqvdhrusefrxqfkygqcvvgmbapqdauvlddltrthtmdgjyuyrsvfjupelqpapqczichyycchmfecwuahrteazovkglwvfsfbvteyspycwqwawtzozmuzloxmqgfzddjwezvwrivbcyopnafdacdwwzvihmuwrrrojeusoeiilovzuymocyeulqaclyqnnblwddnrhwirbdqfytpfenkcbdydwyourosuzasvskgwwmzspxeoaxsdhufyfsmlkpqybrvkslvywmbdwzcuevxwnfvwqbhljmbsuganwzaqatvnooxdkqtminicxrtnxbsaovbvzhklfqhqcgouknbuithffgyjrdiyuagjtfhnzfxefarmcsjblexflvsknllviavdhulnxvfpxhquyyhnhqkhsphprbmwregqodoceoxoctssflaotzrirhhsfykpwwywicwnpriphdphlfgldomlgxjjganiyqaoqxoawpnthcdwbbtwtccdpszjwwrbzntcetxbmomqyyylstirhlyiwopfucwnldpgbeqjpzhdnlsigaozxqaaerfpxujxyireefpqknzdxcvchllcnazzfzwgluywpxyqaqhmuazsunsfzessteyquiqkiwsavdumyyxrpkmahyojjozlgvntsjitcurgfqvzmrifitnijntfaefxukjabbixumyncxkupzxrijdzuelskydnvodwzcqxhrnbeybxkwxinbcrmzvvbgifpaiyuegsijormxtynubxvpglobqznpesxklgyuyibbhpnwecvgqcghfgqhstkpiiinxummrwoxzcbknipwoaxecdvtjaztjbtzlibotwaufnqlvfnwwiufrrvoowrnvclbmuiomphxgpzszlsrlopqworpksuflzrrebncbkbqogvnizpvginfjovyhuhwxjocxmsbtywjxnoobypxxzmqkwajfmqxzsbdkukqncyxnusmfposombjvkcrsqhanfjyaebhueffbdkubyehlwmlctjwlglrbglopxuojkhmrgmdatseadnbmuvzxtriiwvlrppikvjlkizkyijymafewvdtgbtvradyiwdjhskqvwppuocflzkviwkzkarpgpgywlepqfewjyxgyznfpbwswulqptdvmpdwgthyrskpmtdctagjenbkyffpeqdropdzgqokvmsjuernpijnrmzghppuhbqlzugmjkshpnoaynwmbpcyftiiesiwkcsihgangxcimivykpebfjzzwsxfxuknsejdgorwdvrnyjyeocqohwblojsdbpqfeqxjlobnfmvkqglgilyicdmgnkbcldokjtumzxliwasfjoqsroftzgfnzjoqurnlzoexujdeutjpnyzxhqgwczwolrvoxzhnxzbwmfhfybkfjdjihwyoxpfhvzbleyfjcaafrszpfjbznpidkkrvpfedqbtfsjirdcdaxvnhhrfztmbkepjqgoqycojttshtvsutheyikebdqovpkjuntfjkhfaefzfzuiteiavjmqqzocxqyumqjrryrzhvkhbsueavosomvpyleyyiebbpzaarfqjpqnxxmoevrrbwnacyijtqotkwhvrbdflydcumxstqajnuocfjcltqczjqnvjrviyewejexifrqcnebrgbompjwlaiawlnhcixuvksnymvcwaxvanvypxmdmlsnvaumilqcgtzozkzgkxncaomxkqntqntjycpgayyxrxvracnfpzzqgovpwdijgmscdsphclttqiwqohwkkdenforbsgnfjfdlvigdvtlwgekevmpkyiurofedqjyxuzbibcnuwgqzbiavepacgbeaapgzhqwwmbttacnkqudgyiabdkimrarlpzhmaietevmvphggjwmowdryltcllvmcwbtqhlkukuuglonnujdvsugtprvvubkfeerwkqsqyejaltvplmnzgrmybpuouvnuqvibwfwibrkqmmxeqxjqzejpxnaykwgoxavqqjwmaquehcmvpmnhwtvdjqponvfspzblgrepgysrxxzwsgkbttmmytekkqpfhwgnhcyhkjoavxoqekwhtdqvlodctztrytywmiwamgbmofddvbfzongifcyozbajecbhdtzhmydjpvnbgmvdthbidajnlocmrltsnxcvbizdecjgvrprxrrizoayqptmnvlwvmemvjftluzrxcsjhxyheshlrqhbkjjupngvpvnqxnerfvanhjrabbqjrnghyhvrfqjivipffnmplhxsdhzcpdihnoiowdfkwwcuahwpjjqzbqfmsrzixkaawhvvbcawdgbnruhzdxhgatnvsfuqhoxentyezcceojfxqpwczqmslxpfraxbwcblhhbykdrfruqlaqgjimblxrealosefenelzwyovliotfptswlszhyyphmnufthpyskcpfzfesmzgaiqpfjtexofwgowmlwyqlnryabnhsddfzerpvavwmpjdvhsdbzoxztuexenxyaxvedkoznpiiezzuolijsrfthssrxznzyrqukwqvhixfmvbzioyghguthymrhwimfchzpnlmhmegwsppznattlzptqerrxprcqwzpybupjrtxtpykowssxdipnbogyuiqoexomdlrfxixcqhgwpvgaomhjvpdrzvsnobtaioygyvtgxxbvcxwkeslukyunuiptzixqouvcyjqkdarbfzwijvnvfmnmzavuiwkdkljxlrmwchqzwobtyhqbrinugozzuozykicdqtpgylgqjfgwbczedfbnokhotpxuaymqplwkesutffjeufihwtmrbqcfgkvgnkacdcwaljwvhhahgdhseqogsdfbhcnedkyfmjyatrznlpninpwaoplyugwsjjaeqqsfkssyhcslzqivdgngxfzzmcifjmtnjsdsmfejxrjjckpblykticupwvdcaeuaemdngjuteeirjbhlfqqsvtyhnfrrkiczygcivbvbkyrjutaibxwjlwolsxprhbeqcetmrevsdraqmoovfbcggesuxwzuviltkmoikoquxaatzgpsmlbmxwzfvaqcefveiavjuucypfbdqopsvpywshewjqdzahhlkhfitudirgqtvkdwzhxzeinjhshelrnqbkghhhzkqighliovbdzncqqatzugxrexkiaqkuzzpuassitqrgspxvnwwdtnjzwaelnnlnlxoycabfsbswgypxvmllonkcwsqwaovhvzeuxvstgopwpuknqqeakpirzmuoyutjyhqctmsubljrgdpuknbwcwgmtiztucufeanmbrpcwlbtkveozdlgajokzzadwnhqigklfaogbxkrrjdhjalfcnziehwttmxqryxyrqzbibqwroxorbcecjtfpofeloivhrwzwjqrztekvqxddxyjbnueurzzplswtpebssrqhromuhxagyfvjlaijrhxpppzrbqvdodsxhbiogspugvujxfjebqxetnaclgpfhzflrjkkpfojenfimspdfvrhiqjuqvysojpaoxrnjjmtwwdpducyjyvsiwevsuvmbuxgpbashmnshddkgjdviwelqkszpclhirrfsxyxncyfxfcipiemhrsnxrhdfffophkyjoogqawnvtrgvmbjrymzddtvnwobzrrwulambdsnxfdagdgshekhddyxujfhjcxvfrnbgeveeqtqvefrtrstgmtdpbifsrycqopoxlevtpmjvgpertlgcxlsylstnnjuddwcxqclzxvbgiemzmtgaygnencldlarjjzucnbwjfqetbuzqnuwczzbxcddyhdpehexbtyekkihluhpdfhslaqeyjrwkynhzdbwngskwjcjlvmgvvbzqzpulcwltmfdqlmkuzutujnsnrgibeextkunfjeduxlvvsimfytvoekdcsmhgysfigtcswbrfmytwudxxyfkrwfleqcvlkwazjltbppjuxzehzzpmkjuvglarxreptodhmskyvhdydujxtpyhempsupuiksgqpbjuxfcyhehprnytibyqwqogzmiobyhzepbjckgbatbfqbheyqgqdbhghpvmuicuocobwvuuvoyffgkaykhzuqskbadzxguyybbbpelrbzqitsqjpgezbewisndykpfyuuheghzhyvisjdrmzgyectrlnketlmkafhtlmbilufjpnlrybexklwnygphweixyfuerovkwktqfluknzbchpfkpqrbblfblsmcfakjcqarfoavjokjobchxkajqfqxkbdffskxmiauxwzpbomdjtqrdzkydfohexhbexakqkfboxnfehgjvzfhcnpubatmvbztybhzfisoqedggnijhfnmgxeyzzapwhemczfnspphtsvjqcjasyqzawfryxtvmqntwcevncgyzqqfzlnptbgoylxxdlmfbaofopnwiuhnfbvywglkdaxlazyazrbmsbsbokaxyjiyqaawmiymhjhkeubplngaqlxgattqlzonwkcnpmourvnvbizvijahgprkkshpoqgwkrbmialoyqtppltuoadtbxaanunqvzgcpeajqbtfdtfocgecyixtdkllatncmjsfpcwlgtoyhxvdlvogaibcunagyxjvklnmqkwvkmfcbvtybysgqhsmwywvaaldbtvjvlwtrpzoqpibwaojmevwzbjnmnpurssofaofwsbdcvvxoprerbsiquihytmcupnsexhsldfdqigxcjuzdqfrtureptvsxwijpoerehcpvixwntlzchhfvuftnhkvtsdwpioepviekxvekisexeafzqdltsabsclgywjwbewbupdagtjvjanmqgwqglmfawzxoahjuvllanavivecwtwbmccrnhlphcgoouuoogchplhnrccmbwtwcevivanallvujhaoxzwafmlgqwgqmnajvjtgadpubwebwjwyglcsbastldqzfaexesikevxkeivpeoipwdstvkhntfuvfhhczltnwxivpchereopjiwxsvtperutrfqdzujcxgiqdfdlshxesnpucmtyhiuqisbrerpoxvvcdbswfoafossrupnmnjbzwvemjoawbipqozprtwlvjvtbdlaavwywmshqgsybytvbcfmkvwkqmnlkvjxyganucbiagovldvxhyotglwcpfsjmcntallkdtxiycegcoftdftbqjaepcgzvqnunaaxbtdaoutlpptqyolaimbrkwgqophskkrpghajivzibvnvruompnckwnozlqttagxlqagnlpbuekhjhmyimwaaqyijyxakobsbsmbrzayzalxadklgwyvbfnhuiwnpofoabfmldxxlyogbtpnlzfqqzygcnvecwtnqmvtxyrfwazqysajcqjvsthppsnfzcmehwpazzyexgmnfhjinggdeqosifzhbytzbvmtabupnchfzvjghefnxobfkqkaxebhxehofdykzdrqtjdmobpzwxuaimxksffdbkxqfqjakxhcbojkojvaofraqcjkafcmslbflbbrqpkfphcbznkulfqtkwkvoreufyxiewhpgynwlkxebyrlnpjfulibmlthfakmlteknlrtceygzmrdjsivyhzhgehuuyfpkydnsiwebzegpjqstiqzbrlepbbbyyugxzdabksquzhkyakgffyovuuvwbocouciumvphghbdqgqyehbqfbtabgkcjbpezhyboimzgoqwqybitynrphehycfxujbpqgskiupuspmehyptxjudydhvyksmhdotperxralgvujkmpzzhezxujppbtljzawklvcqelfwrkfyxxduwtymfrbwsctgifsyghmscdkeovtyfmisvvlxudejfnuktxeebigrnsnjutuzukmlqdfmtlwclupzqzbvvgmvljcjwksgnwbdzhnykwrjyeqalshfdphulhikkeytbxehepdhyddcxbzzcwunqzubteqfjwbncuzjjraldlcnengyagtmzmeigbvxzlcqxcwddujnntslyslxcgltrepgvjmptvelxopoqcyrsfibpdtmgtsrtrfevqtqeevegbnrfvxcjhfjuxyddhkehsgdgadfxnsdbmaluwrrzbownvtddzmyrjbmvgrtvnwaqgoojykhpofffdhrxnsrhmeipicfxfycnxyxsfrrihlcpzskqlewivdjgkddhsnmhsabpgxubmvusvewisvyjycudpdwwtmjjnrxoapjosyvqujqihrvfdpsmifnejofpkkjrlfzhfpglcantexqbejfxjuvgupsgoibhxsdodvqbrzpppxhrjialjvfygaxhumorhqrssbeptwslpzzrueunbjyxddxqvketzrqjwzwrhviolefopftjcecbroxorwqbibzqryxyrqxmttwheizncflajhdjrrkxbgoaflkgiqhnwdazzkojagldzoevktblwcprbmnaefucutzitmgwcwbnkupdgrjlbusmtcqhyjtuyoumzripkaeqqnkupwpogtsvxuezvhvoawqswcknollmvxpygwsbsfbacyoxlnlnnleawzjntdwwnvxpsgrqtissaupzzukqaikxerxguztaqqcnzdbvoilhgiqkzhhhgkbqnrlehshjniezxhzwdkvtqgridutifhklhhazdqjwehswypvspoqdbfpycuujvaievfecqavfzwxmblmspgztaaxuqokiomktlivuzwxuseggcbfvoomqardsvermtecqebhrpxslowljwxbiatujrykbvbvicgyzcikrrfnhytvsqqflhbjrieetujgndmeaueacdvwpucitkylbpkcjjrxjefmsdsjntmjficmzzfxgngdviqzlschysskfsqqeajjswguylpoawpninplnzrtayjmfykdenchbfdsgoqeshdghahhvwjlawcdcakngvkgfcqbrmtwhifuejfftusekwlpqmyauxptohkonbfdezcbwgfjqglygptqdcikyzouzzogunirbqhytbowzqhcwmrlxjlkdkwiuvazmnmfvnvjiwzfbradkqjycvuoqxiztpiunuykulsekwxcvbxxgtvygyoiatbonsvzrdpvjhmoagvpwghqcxixfrldmoxeoqiuygobnpidxsswokyptxtrjpubypzwqcrpxrreqtpzlttanzppswgemhmlnpzhcfmiwhrmyhtughgyoizbvmfxihvqwkuqryznzxrsshtfrsjilouzzeiipnzokdevxayxnexeutzxozbdshvdjpmwvavprezfddshnbayrnlqywlmwogwfoxetjfpqiagzmsefzfpcksyphtfunmhpyyhzslwstpftoilvoywzlenefesolaerxlbmijgqalqurfrdkybhhlbcwbxarfpxlsmqzcwpqxfjoecczeytnexohqufsvntaghxdzhurnbgdwacbvvhwaakxizrsmfqbzqjjpwhaucwwkfdwoionhidpczhdsxhlpmnffpivijqfrvhyhgnrjqbbarjhnavfrenxqnvpvgnpujjkbhqrlhsehyxhjscxrzultfjvmemvwlvnmtpqyaozirrxrprvgjcedzibvcxnstlrmcolnjadibhtdvmgbnvpjdymhztdhbcejabzoycfignozfbvddfombgmawimwytyrtztcdolvqdthwkeqoxvaojkhychngwhfpqkketymmttbkgswzxxrsygperglbzpsfvnopqjdvtwhnmpvmcheuqamwjqqvaxogwkyanxpjezqjxqexmmqkrbiwfwbivqunvuoupbymrgznmlpvtlajeyqsqkwreefkbuvvrptgusvdjunnolguukuklhqtbwcmvllctlyrdwomwjgghpvmveteiamhzplrarmikdbaiygduqkncattbmwwqhzgpaaebgcapevaibzqgwuncbibzuxyjqdeforuiykpmvekegwltvdgivldfjfngsbrofnedkkwhoqwiqttlchpsdcsmgjidwpvogqzzpfncarvxrxyyagpcyjtnqtnqkxmoacnxkgzkzoztgcqlimuavnslmdmxpyvnavxawcvmynskvuxichnlwaialwjpmobgrbencqrfixejeweyivrjvnqjzcqtlcjfcounjaqtsxmucdylfdbrvhwktoqtjiycanwbrrveomxxnqpjqfraazpbbeiyyelypvmosovaeusbhkvhzryrrjqmuyqxcozqqmjvaietiuzfzfeafhkjftnujkpvoqdbekiyehtusvthsttjocyqogqjpekbmtzfrhhnvxadcdrijsftbqdefpvrkkdipnzbjfpzsrfaacjfyelbzvhfpxoywhijdjfkbyfhfmwbzxnhzxovrlowzcwgqhxzynpjtuedjuxeozlnruqojznfgztforsqojfsawilxzmutjkodlcbkngmdciyliglgqkvmfnboljxqefqpbdsjolbwhoqcoeyjynrvdwrogdjesnkuxfxswzzjfbepkyvimicxgnaghisckwiseiitfycpbmwnyaonphskjmguzlqbhupphgzmrnjipnreujsmvkoqgzdpordqepffykbnejgatcdtmpksryhtgwdpmvdtpqluwswbpfnzygxyjwefqpelwygpgprakzkwivkzlfcouppwvqkshjdwiydarvtbgtdvwefamyjiykzikljvkipprlvwiirtxzvumbndaestadmgrmhkjouxpolgbrlglwjtclmwlheybukdbffeuhbeayjfnahqsrckvjbmosopfmsunxycnqkukdbszxqmfjawkqmzxxpyboonxjwytbsmxcojxwhuhyvojfnigvpzinvgoqbkbcnberrzlfuskprowqpolrslzszpgxhpmoiumblcvnrwoovrrfuiwwnfvlqnfuawtobilztbjtzajtvdcexaowpinkbczxowrmmuxniiipktshqgfhgcqgvcewnphbbiyuyglkxsepnzqbolgpvxbunytxmrojisgeuyiapfigbvvzmrcbnixwkxbyebnrhxqczwdovndyksleuzdjirxzpukxcnymuxibbajkuxfeaftnjintifirmzvqfgructijstnvglzojjoyhamkprxyymudvaswikqiuqyetssezfsnuszaumhqaqyxpwyulgwzfzzancllhcvcxdznkqpfeeriyxjuxpfreaaqxzoagislndhzpjqebgpdlnwcufpowiylhritslyyyqmombxtectnzbrwwjzspdcctwtbbwdchtnpwaoxqoaqyinagjjxglmodlgflhpdhpirpnwciwywwpkyfshhrirztoalfsstcoxoecodoqgerwmbrphpshkqhnhyyuqhxpfvxnluhdvaivllnksvlfxelbjscmrafexfznhftjgauyidrjygffhtiubnkuogcqhqflkhzvbvoasbxntrxcinimtqkdxoonvtaqazwnagusbmjlhbqwvfnwxveuczwdbmwyvlskvrbyqpklmsfyfuhdsxaoexpszmwwgksvsazusoruoywdydbcknefptyfqdbriwhrnddwlbnnqylcaqlueycomyuzvoliieosuejorrrwumhivzwwdcadfanpoycbvirwvzewjddzfgqmxolzumzoztwawqwcypsyetvbfsfvwlgkvozaetrhauwcefmhccyyhcizcqpapqlepujfvsryuyjgdmthtrtlddlvuadqpabmgvvcqgykfqxrfesurhdvqslneyeuwrvhznjuauvchdcjeokxcmgqvdssqycxczsojhlabqbjxtjxivgyunajidlwdnubtpirjhjpefjujwfxabeauvheltgykouqckkqbhiswxbbmsqiypvejanfnzgnlsndgnvwutnknhtleuuuaqqvpiqzzxqskvodjqkryvkajmyetjeavqqozcovtllcitacfkoynzzjmckmchdddbkjwpjgzqcsrcpojzxffvpeomwxwtqmjfmxjmkbwsgimhgnprkxmlhziuoateossqkdcnxvemjmvomqaaqbhowlwgtwbhyzyejejtwyukyufvikkoahabegmushqujpvzicbjnltpihrhaweqvnupqzvlpaipfldzxeytsvizluhqpammrcduswjllhqcgefaigjsuyajrszbzkncuelxqcdhjbgxolvktpoosiudkxveofzghyipijjthmcxnwvejeimevfntqhkcazqizqyxpuvhzjdkxhtnhtlxoykfwhipshtqojgeyvqbkcghbzeoeejdywthrpfkfrhsmtjxubpmuobzrddavuoxlwegawupxfyurachwcnjibvtnbsyvbtsdjyowuvxehldiqtwycopnfxkqitlxhkfsrysshhsiuflqddtqvgigbnnlyoxqefegayuqdkqemcacchqdvwbkobobvnmodiyfgqmajyglkndfqcvxyvbkhwpdyeneaccusuywywgmnvpebahdosyxtemmcuibhdagzaujrupszsfonijypuejccojidoeuhqkcxoeurofibqbrlvmctaruafhjjyiksdlqxlxeybxipixuksszraqewyeodzwbearancreorngcbcmaiwnsklsilvbwclucfujiwuhzvdiregutztknfrcsccrraeoyeivdksvxikrrbvfbususxdfibhxwocyazppasjvicqregrltdtsclnyqqbktnppntlkwmizzrsywklwkusmjvaslnihwtzmglerebvkoziljnhpnpvyyojbouvqvehvqyqgevaecodprreanwrdmauhxpiklwfrsxmmqquxlkcidcstuhaaeqrqpwujbbqycsverxysmchqeffxpqkkvdaarqljdtnnbuwqpkusypdovpkxoxxxcasflhvxivjdncsvbxqtfkrasdvsmddnqefbqmcijknrkifrrbtspgtywqnainnqjippbzgumkycwpbpnbplokakpyqtrbzbbqpdnuxoevctwviqohsoeapsocmxggjdyjqeebsdawtfjmdskmuwimrlljiyqeyicgqlkpglevhogsvqhkyjmgocyzddylohlfecvkcsfjydyxsipxryhiadgzecxshzovioorosjlwkafnvqwgfwhldltwngkxqpngslszdyeqxfayujvaozcfgbuuncswxvsvvqvlukaoyoakybhosgfmxgxbdplyocywdlzxouqoaatsoenpixotkjgqceebytbvxkukivltcounnhooekgtnpthfjfswfufbrjyplvlfjlmxcplkwsiotkgtxcbwrfmvzcmdkzgevvlsnyoaalxosbjldadfztczmdowfhzyrctetdhckispjwxiplpmmvghhuxfxnogoutbtoqokxjywgfsceanaqbgqxykthxgigtgltwrwmdjyqrbcenmppwwdkqezyrhokyeuodexnhmtlbrywcgfkhqbddmhenitsfnqsqxduohngogbspgouquqgycsexwlsnycntnjtxgnobdevybhhuutmurdphtsizpnfbwxmhlquirmntjslkshqmwekfrrjabjinspfuaixahfzlsdcbommfnrjasmjrchjpbmregytvbtxcqwwccccwonoodnqmsqeoscvmqkpetahrsuenlypupcofpgdwoqsoofjctuhaggpybkozqwjozokcwpwebmmflnrplcrjzeyvfmfuejakxuuisqksbcwdnoovkbbuofzlvtywntjtkytnuedywvtjoeenkfvnyvrlaxgfneeajojpcddetwoengxyzilxhshgykneuelkxavrlgxfdotkyamjqlcoxdqxcjcxchfomutjdayxyjicshwncnrzfpdtkoklbveolamhujrlovqeedkdvqqrdabgowzquqbqxbganszrtfmtejngfbcbjzaxblfggqklqoeknzgnwhpcsexbmrbtgvpsfwnuhlkqivrrozckbfpexpirtjwxutgtdmrjwrtxvsyvfxyewvitvnntkjixegewgxqahhtsnntmkexkzhogrvgoophudduycfsljnmebayjbrvptpuiganevgazvcktmtfbzqizkfxiavxmbslrsjokvwkzjwqpawpzwigvzqryiguwqyuyajfatfohfxhxvqhjypfyqjspdpkpvemhlzdovmouzaffdwmjwfgzwtbmrwripjywwcdxznapxtuxsltmswosuogdetihopaefsksruosqoeqbffyycqbygcvmosblwnspskagleiduywbvajahzlgfntfodqovgnpqzfjqauxlltcffnbzsdupeklhsugzqyvicmrkbiofnimggeqhsuxsgopenrjvbirzgvhjgsofvrkqotohqcfbxycmeeriabwcixsmeiccsfdwnozaihtzkyoeuggjvvurpdabkssqtktdhhdpmekffzvanlrrjsvypimfiwqpxvhppoetqfpzctmyvxwytsyhxnejvwamcrprxyweqmxzmtznqwjwqsobvgkqglrkuicinrlmuothemllzsjkrkhqzvuhgjsbihktquztlgrsodibfdhwvfdmswibibzbcbjeshcnsceqyjcpuphplqovhvfbchpxiwhhgkfftnsfgvyhajfhsyieegoqnsvutvisupjuzpxvstljssqaojjmdqgdzxwsxmhwbhzqxmookhfuahggxnqrvhqglbrommcjwrvhcjpgifljvghmdwbghikrwqtmvxkjfoxnngoclbefrchwgjfgwepgurybwuguxreeqqodfsneokuystoyghootophuxudiuajfkviodovtjtslihpduqwgogjhsgrcqwsvqgnzkbrfnayleqehfvgjckyvebmpstcuheghlhnxzxeldlldjriigobpmzycyjppopyrxsmxjxxsgeshzzrhdlbqtsjefdnarlywevadbywzzbdsrgvwvaovyjnyvxybrjxcrbvfnxrfqdmhpnsvkvebhofmimbylwakpsygxyyzcvegdudtaavnpftktqkifavsobsmbfznrzgtmubyhkwxskotzdvaqeehcmtnzvdyfqvtlkuoqwnmxalsniherwkyndylmjxipucnarkbindymazmokvdhpnogcotiexxkdkpgoqndoswcjbyvumypwrtwohuxmgpmtmufjuwioymgzplcybjtypgztsbtcxbtosdsgfbacgsftcjpdzxgqissvipbwtumkzrvisirtabspbilbnctpqtlibkvhdbqvxtzguhlybvonerpttmfztjfcdcebqbxdpuwyuipwotbgjgxbxuarinnzvqbksuyncyutuocdzfyvqryimsbclzmxqtomuirkyqgrfqbkypabhmtvuawbnytlbrgiiadzhqydlkjrtzzpomivdizkcccndqmfhftojkvaqsktetthndhxihbleeoqtijncypdculgbrgxrbhkyptdurvwnsmbdfqbbbvjnmgtprxgwplvbgflszvaotfxemsbqvusnghuxdxabkdtpgmbbaaachzbkinaeokpzqfhmshftpxezwwicsazfqrvqxygnoiqqmalvvxkhsvgrqmwuqoozjzwpkfwehooyzndmvnwdekbgivtpvaeiacsdxhosnciuydsxkojyfvcjaoowlimahcxfcatxrtzuosghfeybatlwktqwzkvdgngngolftbrtkwmjrytfipnvlzjbrfgkhjoptoiamdezryidpxygkcfeykouseuiaunnrcswfstahgngmxihlpertgbhrvteacpqcgxxsresbqjiwnudkrmwrytrjbvcnjaujdnutskcpcjqulorqklecmjyqyakmgfksvkeafazgrcmwqxvvjxbqputkvztopuptwptgclqehytxmduxcjvdbuyikwrjmbpnblekxtwxoazcxhhzwaitygxrasoaicqjaprtcwuyghpkjgcyiyggliwtiernfhmvzymalaleywzfffhafyewmbzucneqntggheohlnanugqhzyjedwbbgqfxuubchyvqirjaowsffmgzhhrzugecukjtjsiinphfdxnbgoavjdllmobigyyptggzxdkahxrlotwoicinvqqafecipncrkwedovjiugrcmulqikrzhipnbauiijhimdigbemwibchevvanfesvgjqfdnqshkpqsvbbealwxukebzfsnbtlznlbmawpswturfrlsruzgytmccxfzupzfstmymcdklmyyfqfyrcncbofbktoetnojljeicdqujhwqyqvuyeiobddhktjdbwtgjbfkijqovlcdcgstovqkfyrmtesrxfiythwsietwuhceqyuhwusktzjdjxcbtxdbtksnipjmzirajsnivddvdpttablreqjprlxevpgghnizvqmwmxfrcecxnwpzzhejcnwlwkoyvkiemmtvuizolummjhlgidmnzywwdhzpyozyqvcvgzqnryripcrgbogaluzxudnxvezvvaaunaxfmsfgooimxvqhratgqgqdzjgkcjkenpvtdvhuseataerdlvupmzuzawerkzsncpzptixynmlfstassizrwwgewgwebpnggxwudjfkzhdofdnkwhldaqswkqandnqmimzwjogjgmhqpsxepinmfywunaepdriipstogvjzbnkhxxcwvxsfnbdblxrduzgcztxdttajainmibjbqploxdemlydeukzxcmqnobpywqkmljklwgbnitpqbomcgwygyhdjyxrbicugqixvszazogvsbzijpjslhtytejafkpdbzhkkfkxqrylfgzwlwknqgkfggiddekyfkuoauqiswfvtdkohydbvvpskjvffgoqexaancafbekpanvhxrgeskdiksygkbgoczdcxmguorzydozwdrihmeumnbbetxdjxbchgddpmmzshcmnoglxquizdydknfzblnazocrwxwxcnqrxrdbxzsmuhokvhnvoiwbgwlnebpfadxdjqgsosutxsojwkhmtjmdxbbugxogyegebdbsyktnjjikytgxuplckqyjqidesucmqzhhzmwcyegasjiulrqjwezadezfkevpmernyoqecwedabmiduqetcpnamjjszigqrbypcbzjmszuajtqeezbcytndftirfpdbmkdumkawmgzmquxqgsdeflvzmslzjnotvfsvmrvfyrnoiktvaxckhixshmclpewlgaxqfszpdmktonmgefmstuycjciyhvihtyokfmwqcpfuoruwebubikjlkjxuwjagnoqopzicbfdmzmbftkhvumflzmtjcwcrsmtprpoizwickwpeqfmsjmaguqgntxsgubaeqaxrobpazjrechnrjqabizmqjoxxhhmfutuwusdikyixcwivzlrfucrjtonkurhqnxdejfwrqbuatepohnpmswgpdwzkiirjcmodxpjjwmbfxbpqvvbkvurvbfhabvszwhiqagritjxaczrwkvyhuctoysnzqntmhdaslpewwqaygabszctvfecdfbcvobvcezvsoqvsxdjhdbelbqmyqqorrobihbpdtohxkvvozsvmsivqioheznkjaanwsletmpfdwwxxexozcsaixpuvzruzlilcmwdmzxcywkmthpkjryxwdqnruiovajcytqphynipnrrllthrtbpxdklupcsjnlylmhzvvagdszvhmcvksasrxucrghxorcshqlguqrwzrwjqnlcfzripclspeicomkuykhptnuqnbtxrtlkfshtyrmbrfhveecyoavaanyjcfnsuttkqxuujafglqfeyviyqbtzvukphaqsbfjhttvmirihuekkdkjyhsciuuazjzljuaiovwtxhhoxjhkagwamooggnbbkxpwebgdwhfhydxyszpgrepqpdkomakazfxcsmvbzkoqjmbooeagjbrucmvfnwizfcpsuyrbnwgvutgppazjvklxubjvntgnulbxbslyrjpizygazjdyadhsrryugphjoeegjrqlohkkhitcnyzvokvafgveoamsodczhjjvugxlvmwiuojtpyojmkfkyioehzpuvfvmhllcczilhqqkckoisfotbwurnpljayaziiqqjvghseknoadpeeaxapxpqpkzdyrprtrvgkrmurjeskzgsvkcexpmxwvawnsmmyvhqzslwmljskggslmoyheuicecywqrrgdejzcgwbzwstxypbfyxlefejxdfbidyqgfazsebflmsghdwuvtaporsggdnpokeawjmwlwsmilnlmrnqgrdzpllkehnunqiyawcwbhvequikvtmngicllakrhdrwjrxudstlfwdalsgxjudrqvsgrzjuxounjhlsbwbukqiiaqsryhtqysaielmfuqqnenysaumqvwmqiurtelzfwruvbvmgyexcpdzqwelvvoanukjaglifndebziepeerrbnuapjmerviqdwrbcnaxplbpqycxvbrirwsaphwumnscofkxiffpinrzeyhbyiucihrwlzfmghhvzpljtrwdcipzobltyvttvsfevgduljgqyqnnvfzhofkakbbjofwxwdfnvrccbemfudepbcozunhutgsmpquabfagryuatabiesfufvczexgzbgqbgjslxqnfglqoevdmlvjfcosqwvepbbmncdbxxtqjvzfubhaitzhaptnhxdmgrrgyubalwircnuwhhouhddidtxfudcwcyjjefbkndticcwntjfnlsjabmddtpcmzrwydponwinrqejgrhitzeljxcggsabngmdoowqcmygohbkbwpubqpqjdqianookiskxtmegfkotocnclzaznwohqcoeiqqorjaolhhmuppmhmsqxgdzicocuzsqotmbpfuevcrhseukuopljcghkpggrjinunyurzryzpijdpgezusjtllhdbnlywhjlttddmthjzovrahadrxcceazzkqovquaeuytrdwmkwpljxclszwmuzutntclbgwrizamdyejqxaadkfzjinqvrckcyhmmkghtpzmvehstiwswscohgbfxtuthebgrujfhbhxfwfswcchwdvuegmrobcmigfzfubkgwoaiopizngcpciambrlwraakuntzblzumasrygawxaswqxrefhlfeswptokdbcafqaunikgpuaoxdlbhfzvgwugykbsjprzryyhjvzmviyywoyalwzyruscxuwhfcdkrxbeexsmxgzadbjwouhunrijypsvrzeypjfmovpiczdhopjuggzklquudkxubyvxcqaegjbqlrxusgycfoyfkyqcuescpwahmapjmmitlaireuhsopwtlrkrgxdzrmlodkapylbzevgypsfhothdqtajqrsquyuqftjhixelntobgjvovmxkgqzunjnptyyjgvildzkanfwpvvihxjgszctjekjscimvxhndriayfdynhfdvqmkqmnanvteorvjwuqqopfubbkjktwvtegzfkpgpifhpmwsfxhwmsqlkdccmctvqosximynjtyjfrxfshpgmzuuqymgymqbehfgqimvclwngmsheqsgmininrpkgbgidkbyfjgywyidhycutmniryhsxhtecvlnhfjqfaebybpeqjgqmvqucwobrzescsnglrcxigdjuanluufdertscfrdrlkqhjwanhzdfkmfqobjzlpzyaokzdevaibrxfpodeifrreloymqypgmywgiowniooxdfypnsfrgqlqekgbzdtesmqfxyckvympffudvyahlzpqgscgkjhokpdwvjdsaddhijqohpzzrrbjmoiwpilobabgykzbvrbvbpieeugxdgzdlisgpseojbwogxnmsbdwboffcmhexpvkkdfdyoqywhoxuodxupnbnigcpbmdjizhdbtectmxzjhkezeaswuuhwuydntntrviugxvtyydkezdrbbjojuahxvduimrhfisetcebseszpqidnraauskkfywpxneizgoeckcuhvbxkhqjznawqowjtctvxjdhrcqnqbvitsjggnhcmmclrgwuaduyqxayrhlyjhycthkohnklxovsemjqnusrygnxelqdujsxhmprxpoaxxtdnrrtudijrsutgxccgzxqnttghybzbprdjzoulqajwasrpqdqcsscmqbfcimsedxkyczxoitsaxhblcdpuecwzhtedgsrulcsndngcuijvxinfrhenfdgynqiykokpdsfsnfpydwkgseairskppdbavqrmoukzoussnsnhlpzpldvzeaaoygujiqrnkmgvowjbgkzsalenjfxluiqjtrwwtfpdnlomaiumzwwxwxkzpykaqptonledmcuqxdgrgegeetlemxqltesfcqtvpiduyslgoqqhmyulgwzxqdaxwtayqkmkxkdotlxokkwlqvudlvenoucwmnuxgdehcleltwmaqawdokhaywdcqlbhnpxzchhkwhrsmblicznamjgxuhsqzqymqlbdesnvvhxzfekseoaqteactcbwoxosvehhcrrsvdykibrhjjzxrrpnslifzpwvbopqghidujhgafmwhuofrcoljzdlsgcmwkeuaeroheioismbhtmsyjvafqcktzmirkjeqdnyilbtfuakjkrfgbllpyzrhzalgbjnwhcbvoysxkmduxtldzkkzdasynrrtvwfmftpybwjrumzsitliqcxjhiqxdlfhsbarjfhgzwpsnnmblcqxuqpsrplkxzrsnwxijcahnkrvcnwfyyhltwuwixvrchfuzdzqcpwjapguoadftkojmrlggnctlrxkquxessakzayafxhvemtkkmifsgleslcarrnminhezubkonrdysnktfuzqsvksicbtbrhhznajtowpuqgphjvyqhuitudxjjpupfcbvttxkzqmxxnkbirnogrlubdhurhjeausuqdsgiplttrwfnpniktvracgrjytnouhoiqvqlbkdggdrfjjgcntoldvdwmnpohvffrhcdwbqhyhpdylpbzfkhlgfpwjzzfdolszmkqcutddqublfqdrcpxxgnihookzoknhoklazyucurukwsepkbszlqlbzwygeqpzjdgqatpaqxqpexraeejcnunlzdspfkbyyvqajyfphqgldtpioxwkfgvrsitusrrhvgyycqwgqnvbnhcqhqgqbycawucvwepjtqspqtitfkwsytnwmalbxgekiznkbxwyldurobvnnogzndkhnxiclsmxtlohckkppkcbdfgskztihhzahsfpbrkbnnviwivyivinxlgmpifrjjjlvmmdzljtkoagckbnohalakkxuhmasqhvrszkzkhsgnuinkpzokkienxngvzvzsdqpswmhgmwrlaghzitqjufpqdulxdbyjeyogdpvmcdcadiwvkjezchhpmwzkhmqzjtcchemlycihuhfngdhtrvsvozxcprtxjcnxpnudbrkinzxxqawjjrkniomeonhaevweestsiexhohmlqmlvcnlsqwyybyqkkkskkgvurlqqtwlgstfdulhjzxolvpaztkoyswkdishwqxmptzejkznnmocpqhriqvncrdmatooaxelingluyexzxgxfvmgopdcnmgisncmgqkpqvspkrpmfvfwmqtpexfnpsrdkruvsdvgzybrgahtsidsksmpocgtuqtexvtbbamegtgtagrtxosrjgulftmhpngtsbywbmtsublgyfablrnjtxwxinawatwyjtckwwijiekaxypfksrujbashnsalgeycqsbfvdjltqehhiqafdeiscldrpzpcdqinrgzmvjxgcnblxyyumkmdeykjbiqlapcngjgckxjcqxlonwvhnyebcmxcaaadtrwfnxgvyfhjqesrutydkbelanqbgsfrltesqoyzjwvmakxmrcfmtyngbrrttnipcdhvvgfyzzzbxxgqgnxkquruwfpemzhkqcabccawgpldtgzgtwjlgnodjikuaotjxbxvfnkjvogbqmrtbkraxxbdylqvxdugqkotjbktbjilhccokajrznmoitgurcsszkambwmmyaoluigyncoflkckpaoanullvjhglpibmsjbpxnjumyeydkwzkzqpeimnffkgtvmjiahuwaroscouqapcfraucdcffdfbqdfizxrgqvipwviutdrdtxwkxcjumroukqovxnuqmfnufryrdzauicazslrfqwgwasomikpidlcvewljlkkbizwndmwsukqvvnrgrwpachfndalzudkrreaqxvykvpvsixexxvweqdrzxmbveexmzpridleqozycudmjmlcnpfuoioopuysrntofwkvqhzcrevfmaumwgiiazvaryddpyipwvinzijcvuxulhozwbmddrtfrmcttfxwdspuctixmzehnrbwolmgtbvdoqwqpfofhmawulumupprcnhavxegofjenqdqyydqjcerrqajigycohkzzsbkbhwsmaukqlhtqvsetubclvmzxksnxhgyxrhdjvfcdbqiveykuiosneuausrmlwgfgwmzvfznosvjuzsyogacstwjuizomifxorqgwpjhpinrictxribcaawchdqtizdtgpbqmqxkenfaxrfzwedwqlucfqhxugfubxkwvaaixlybitmwfyaeipcgycgpikxttcfwdbcyeltpuhdmitneolnbxlwjelhgoppgkqfbinfahaphaiuwnqybkeplyevphrdcxwgbwflvesannlsuejqsckgpprauozpcflbpyzpvgqxtpigxgvcunqlsckokjoycqsxvxxdqrptlbdgnrwxjzdyzxazfduiyfeskyufkmvqgziuztmruddghbitabhenugzmyluksysxizxamkceqickzffmvqgqxwosoclpolunevybdwknhazqyfxvnkfbrpwpgjzixgelcuazmtzcoozmqezmkjsrtocpdraopghpxxdoazhvbodniafsgbimqtkbajgszyaykyhzlthdninikrxuxpgceugbddtjwqwrqfzqcnniiwttnqqjotiqnpcqqvtecaeetspmtsodpjeqhjkmrehtjljpxwrrpdhpwtpxiegartipezxnvyynvecxdwgoghewpmmmrhkvjroojcrhudsripiyptiglqdbneukpkmtnluoaacwkbjkiurpjuzirtcslxkyshcmrcjfjlghvhkqzhtjevpmqaneowtrhtsqeamdmplovzvuolgmifczyetrlaroedpayjyeduyopbxaxouyewlfnvngoncdgtrxwmvswlpvhihkzoujagzfhhlycpcczqsltbagvggcxgtcklqwcdpbyvwfbnvxddurogepzvayupyhdzgyndkfjcqvbbwwfdzmbpiovfguhqqyxczegsrpdkysnzatomtwlgusmdhunnqdidndjwwleqoobpujhapxonmkuyctdnfwhkqejgfhbxhyuguesygtzalzipxuexnbdfvnlwsookdsxngwvqcczhoitbxhsrcetnavakpdpdodfxhfljkewjgomkbhvmtfpplmgslcwxpnjnvkpusyqnorabqwoszpeatgllazkrjdfvyuxtglqanmfpqvmqmxhtexsnsshkhhlghuozztxefrunprcvqdtkliwxbpgwiuatxzfqmifaluibpazwgqybjcogxfzgpmotepvjkrjwdhnryhohoaniwgnaphejdwqwwphwcaqrsdpqtxxnszrqgkzgiomdlfdvaoehgejybkrmnmbkcdppcrhydjdelfobjyzzyvihxwmwwofoyljsmjkkxmkyekcqqaarsevenutvsfezqxofgqzponedlhlygooaghfnojrttqypuipsbudrwgeecgldmgsdnnodkdxpgcajqdfuwccpywdtjanfzwmimizasvoyekguntmavsyeuvqgrtiwhwhlfbyxghiktzrhbgljzzzqzuruwggasyyiupfkclpdkyemgsksurcuskwgjgomzeuwfvbpemcidsiebmrazuuqqvdvoqbbsmqkdgyqmgekveduetnmfqpmqysncmazuiqmsykdhbqfxjygfftpacxnggkkkwkayrbdctyjuhpqqmbkiorbsalqkkzkhdcwmrdpzvmyalrrbgcjslanfmkceqdyljzxxymzbhqaupsmfpzamyjxyylekotpdxyskqmoopcbcqnljvqlgjgiqictdrzuleseuvyrnimtdfxvrirhqdfdmpwwjbzcqzsvlzxhpcoejdynsnyxcrkktiaoplivknnoyihvrtewmznararugqrlmatwgdvucmuwbaugysqrijhakncxfwxcckaffbvjwfkcijabvcogdcwupoggnsxpkqfdvkehqdeoavlvgukktgrmkfnzmwteeubspxscsavbswwcmjvhppdplunrrzcujqgjwmdclahkoaibegjyynosmnhctyxgshheviaecglqyblpcpssxfiuybuhnwjetznkoatnmhaemmvisuubfnclavjffbuhoftdvsfrxqvnicyfyiigdeyqcixugbtinbwtlgteijzaqcolvragkdvpnrknqopzwmloxjdcqqfjumykiqhznbhwwkwlyhhcavhsymdvikbsmotktnuzqspcpxhqqhanqqdctokzmapcewcjstchnxqskvhjoanmdywereirkjasksomklzmytmfijlmfjzwgdnlkhzstxpsflppxeddpsecxkwupkcktzedhpmrtqmmomrfftpamtgxbaurqsdoqrotzknjqwedbtgzigvyqpaybzuozmxuusiybfsuxjuvajcdodnhnwcjnoypouurxfiferzdcepkiwlrlqzuwcdugbnjbslpopoasgmsoukzxpciwxhbmznpsyeuipaeszhgdbfwqphlgruicvyqkrjqoaiaslojxaizkhvwxswowfroouyawjsawdyryfalbfvimhziruntjditmyvudvdypdmhydbeyyqsagwfmgjbdemahluqcrzilweqojrxyhksfmadcfxrhnyhylglerjkdhxvscujttdihulopieflbexczfznwmmcsbcpscdhzutzrywdwdmjqwhzxrsixidcsrwrtrsfjocykrvggmrtrqdbqrtaoglvljnikvmvtgliurigaxzkwjdqmlvzwclcphscsnfvmcfrcazvukqeiysraenaqctiaimkljjinuaxhzarqplsuspqhzvriqioryvabtyrgbeoxqlxnnlebicrmvwfkqznjazpvdzdtervcqqujmswqzskdcuvkbsxzsbfztvnjbiejcwprtxngibbinkhxkgqhrlacbhtilxobntysjutujnjoqvaxvvibgpyyvkilcmuqgddheknrggwxumlnswvdcetqksnicopeisxrmkeaiupebjfphaphdqtevnfzexfxkmygbacdxliujbpholsviyrwjfxiyjboippbioqhkpgxmwqosdzujwsuphoekjzyytzeyohmeongvdmrocmpybakrixilbrfrmnzduauuuvhdauhfpboroexngumdxcrsngbexzeqaylgqrvaiddhchaetnsxqoutthplgpzathrbagbxynurkbfkocspioahbieomhbjxtyqhtmsgdseksrcmzddaycehrdasdhlsvbqrmgyjqwkyccgdjbiixufkbicoimtmscxwkxugnkxxvqmzxqjrxarsaltiqfwpmtmuxzweofnoknvskybgfrcuvbkzqlsfssqxqdktjuakbiwhnlbekinxaoezbtyoufsnlxnrepazfsgrwvtqtpbqiabqrdfefqlstclabfmoesedsbyxhtowndsoxaydpntkxdtzzuzzrfzphtqdzyykmivnofebggkylwqhjyjmocjhahkppufnjfpymdcopftlfmgwfnsknbxebwifbjmdwzcfandkqchlvqmawvjmcigjzefnfluydmfjczsldbuqrqrtnhncvviwbcirwvzuholmivsrarxritxzgmipevfhhvrmzxqrqcxljexzjnhgmkuefufuajmykcwxywacxgdyoczzudzyldjmfenhfywvgbvmmhtjrpjscsfbklzvzbeqdmmfnhpqbmquozatbmgvkpavarafzvsfmgqxbrlclfsqglchllkmyilpusxcpjexetpweblglnsboefuzhmmwgprabhacpaqgaemtfxazfrofnuhfaxecamdwsifyaxfgxwtffykmfcwcsptuecfsficrifocpgvvcsvchzputuvrrpqoyqddcrqbknqpcvtqsjnvvthvqlmbncyifzljrupmvhbsznheicjfgpwfhpsjuxaagjhscfiabyyqxkuevxrftrnilzyguxcopdmrooevryytgtctdrukadpaeecnwgqubfamuyrkcpqxkxmfevkhwpkmqizsixddfixgyvbafkcfblprikiquqlachvtauzstocicviqehnqhxnvduhqqnrxgititlkrotvclukjikwwcpjacfpedkknpzmcatupyncgxmomipimpbrsafwncpzaljvmahvlxztbgxajawfgiurjffinuctjxshrrbkgkdvgajinxepqonhzydywnvmekvmbjhrykvmriwpxxdpkcjxdxmouyywgcoryymidlbloywzwnvoewhecalgqghidkckkoorofhhzagwzaobbissjjlgirgdfzwrbhvxjrlymuidmjwgiruifwznlewnyapepszewtbbkoemafkjsysmtqdlkczseyfxtleaggxmpxbxhbauwadsptriwgkpptzpbbjmwyqbicwlmyukgvomkplvsuhhlttqdvboecabzwgoadipxijkjqbngchzqplyivdslhhtbscgltnemswivflchgzvuhacdisnliodlavvzaswrceevtfmqyhopkjdkhqkqmjykmmuzfgqwwtkvasiqhskmsxgdcmjcjbbkrblskyeekzkajsqhyypzhwgwhddsfckkpqcbrpkfggmbssssoumbcrnpcsyghqyfivmocjhiunsrljcqhbjzxgadfhnnjexvolpspzcgitszyghzncacxfnhzyhuqbagmmehtuvtltjifmcewbcsqpdenwnciogyohcddsjlpmrprvqrqdeutnpwtbnpimeeeohjksrtrzifscztwsuzdgxnwayeuxataocpiqqkwnfpftddwprhxbgartsmsgpmnngvyewutsnfgjooafquapjokoqtglkxsotzslozuwpvcdpmyjwxvsfdebqdbjolkkttyfvvkwpdvlmihhynpehpceteeabzezbxuchflhkehosgiytxwhmtciasfihyxvbdnkovpjrukwajmeepwokwxvokitvuqxrgtekuorgecudhtlosfopzgzmbafzewejxykhvmkrvhsikmobjdfkgrhstboksfbsuilaaxjojxztqqkjkblgicoydlvmqrtrcibrhfbclseigvkycuzkdjfxadebgpftpsyyqjlzjhuxntcfkacfdrnighnlqkxnrmqbhnfcceztdvzeappuupkwojttrzecqybghjmtazijabazdfgaufsaoxqhsoebbligcxfaztkcnqdirtcygfnkpeuykjdigaibdbgfvanpdqnbodybbqkipcluaaozzvktjfxfhubfpiazniqurgtmkhtqvpjhcnsrnorlrecxkvgqlueavxubwbeqkukabximmksiauhdsjwolbctdujpycbkdhgclskmcgflhtgfvzwwzspdktjescraorstdeerdlrcxeftewbttroflhusvqujrlsdigcbrsvxvtlnpjuerqmzocwqzvewxkglsifldalarzfpbcmyfyzhyqumltpttcnolrhbavtwyjjjsdkwohgctbxvwpdmmvrsuhqbaytnyypamyrdgupnzpoyqzostsxvurjhevtxgtuclperzkfboaxyiuinwuvzflwatsdgpkcmgjtmcqpzyzhkafcwtcweiwmchwqjyicogqxqkwvasxkkqzzicrlthbtcsxbbduqxkprilyisunzqbowvuerxmwfnvhjylthmyuqmoqtmaijdczbzazaxfacvmmozmqaiqoaflkebfbfuciyxqnsrnwnsxvyhiakmcsmxsjsofaflpmhytchaatspuhwwshexvcwbwklnszuoshnnpphmmljgavrkdrvdwtielcbkvubetuwavqlshpyhslwswbyqbnhjwytamhenxnffwkwysushrqalmjaoppjqjyfgiittqnxqkzkthvmhicsuhwqxxihfrcvmzhavksgskvmgrndpqmrdvojbjcswbhltntykcdqjjlzzjogbdmigyvtwmvfkqxkkkuwjjaypzphrfgzveyirdmmfxpmbviqcfwywsnvnjdrdjpfjpqgpaouxxfczrggdrtpuiwtayiwwebhyefcklqinkodstwcaxssdommyvzldbteblsxaugcmwxczegbdodfzkjxvfjjdyuyfzexgciratmtryzjaobpjcwluyftzesbqcxqqxpovyaemlmixydroundfruypkacpnokhxlklalirotiflrogmbudnvoljztdejsdlxqfadlircofremnspxnxvqomodebfkqoqomlzzcpanuqnyahxmiqrfqowiaemcgshspuwrrtjwnhegqtbpjlswswgavnkaopnltaeyfcwncqutgrervfpuqfczsvmqlevmhhagkulohwjafevhcgupqbtjvlmrmcwfuzpaswybgjmkwakpmgzigyhjbgksmeemiueirctvdngywtlnzwcthkthwzjnrabzchzjyxtexsjyuujjipfizugtfcoviuktxizjyzydcxzvhrslidglrnxgnrolblronmgszyhlmzkfnuzkqcqkchihrekgsxkrrnzopxnczsobskdkxwcrvzjzjwgbndzsslbdxjtmrrkvybsqcrruvrgiqbecscyqwpzsgrbfdqrfhoghvzwqqrufyaunijoirlhsaoesljdrhugsdckyanjrydymkrtqbuhmfoodbrlornzbunfidrivfuwhzpdfcofwlorgyrllvifmrsyucnjjmzuieufibtsbhbxmaqjgwjkvimaoarrvtxydjgrlhikylmryplhzpgqycxkdkhjzuuapnjtyfvwifxhtjazajsjdbuijsrkzhmtkavyqkvxtrjtwhljjfussxuwvylzwqmcrwgznjczxezluoaymirflplkqxvzzfyiftiuzdwupqttncsexeazhmaplhtfxkwflhprdvtvdnbtpnwaxwkxpiusowcbvkorwqvxzqghitywyifkijzpgwizeojeuwktlytukquqfyrcrplfgfqlgdsjatnhbnejgqefsckjhqkwneguifpjrwwburxqrnnmppzqvsfltagulxupcqjzdewgrthudblnlhgrhtbcmrvrgyecgiksdsjzsivmifljfljcomvsiwiymujuzepxboapajvauyzioysbghdxtpysvelqaswxbguafyrunclonfnvewbmpmutjzkhctolbtpcnglsynrzyjpezxmbbaddtqnwanvyovdzsmmcszfjuhtxkrcuusirhspgoazxqzqdgmcwznrltvewqzzatsjaqsvhkelivxgvaduivtprrkecuppcbyxhojqpnmtxlpcnqsudmxwevkfcbnrcofnxuirynzibbuxhjehgcunlkysoqzxnzrxioyasdmpsaashhkacahdmdkifvdfcckohljsvmxwhipumflxkalpaobgjmatnpjewfeffeznsyrzizkmkjtaialnbmljfxzzsmvgyrifdwqkwaqxmathfbnqntuzybviaisawlkzropdhjpnnlgowoyksliawgscovugkznvlilnksybvtepzbkrtytqkkffclriidxesgddipogcnfqhwpuuvixtxqyqnqngetivnbgvdhkyjnjzekvdxhbnjzzxylchjujtupcocgwizglzfivromdqlarninuikajwvxlqfdmicjyuhmcmngsgwsbdepnjseggpiazcbsjhxscdfvfyunxstkztchrhqvwhczqzgygezrdibnbimxibpbtetowvcpxsktnkzfspalhonnglroomnrebswczvdhjdhpoxvefxnfltdpucjngguefbdtdeydmfhffnmgivbieivlmcxwxjkfocrlcssqjonnsjqlmzocmnezxuubrclymckaciigtpcitpzlkakvmakpsxkhtpctdxblozgwwrercgmkbkykgtrrojsftfnkqkadyvamajiyspullrtabvcmjnphhqxnvhfxeogedkafxvrhirpnhbjcibivgdzdcngwfvrtageshqpmxnsjlxsaogdfzqnbvsrxnnunaogplzaeauucdslcpmbuorfclnhgzhdcgsfqjtdktpduqyuyhzofczabcjowyiqywyyshrlcwjtkdosjrvsoyzyfubtbmxtfezjapxohscpcocuzovhvikjxmrjyixiekzyrbrfgkbutactlviplusqlfabnfahsxokcoujbwndpicfbbbwomikbnpkpbazbbiwhbvqgbzxrpgiakkwufhgyndqeckcnpltbymekiccbztxuvqiieiwcgffewnukukvrkkvpqwcalebviwvcisdacniibzrrncpmaanbhsbqfrlhvwmwnjxamlwcqsaacuppnmcnjaqqgczqwwrlbffxphxfatozrmjomgixkxxgjkrvqnnruzukxlatudauzmgzkitemizffctfsrbvvjvbqxidoxryemkpzmydmpgofcgcttaxxvqftgczfiprsacwosxbwztbndwrzqiofttpvronlzwycvgspukdvahyhlkitouiybmicipzfwhnuonbtmbaxoxqrownsyulcqubfojkuamybbwiiazpkfxnlgsmloncxogcsmuaaiixiriiketbuxmzaphgpqyircvcqkelsqhzdeifudyycyhuzcvgnjjyseoameztfcykidreootxofbiatazlbakqslwsxipervnelwqyjbzsefmsfdibnffisackjgauenxtiqmujsfmwifjwijoaxnbfmwwazrxobrfwckbkalbuvljafpfyxztzhjuzhovekgydwjfgyvoeyfxglkwmmzfclcguopzharhbukohxmtbsvmkiwacgprwgjbszviuflznyejieqgnfbzkwisyroevne";
    cout << validPalindrome::validPalindrome(s) << endl;
    s = "cbbcc";
    cout << validPalindrome::validPalindrome(s) << endl;
    s = "aba";
    cout << validPalindrome::validPalindrome(s) << endl;
    s = "abca";
    cout << validPalindrome::validPalindrome(s) << endl;
    s = "abs";
    cout << validPalindrome::validPalindrome(s) << endl;
}

int main() {
    validPalindrome_test();
    {
    //judgePoint24_test();

        //checkValidString_test();

        //MapSum_test();

        //findLengthOfLCIS_test();

        //findNumberOfLIS_test();

        //flipLights_test();

        //findSecondMinimumValue_test();

        //maximumSwap_test();

        //trimBST_test();

        //findKthNumber668_test();

        //constructArray_test();

        //checkPossibility_test();

        //strangePrinter_test();

        //widthOfBinaryTree_test();

        //imageSmoother_test();

        //isPossible_test();

        //printTree_test();

        //findTarget_test();

        //findDuplicateSubtrees_test();

        //minSteps_test();

        //replaceWords_test();

        //countSubstrings_test();

        //findLongestChain_test();

        //findErrorNums_test();

        //averageOfLevels_test();

        //judgeSquareSum_test();

        //kInversePairs_test();

        //maximumProduct_test();

        //addOneRow_test();

        //leastInterval_test();

        //mergeTrees_test();

        //triangleNumber_test();

        //findLHS_test();

        //fractionAddition_test();

        //postorder_test();

        //preorder_test();

        //outerTrees_test();

        //minDistance_test();

        //findUnsortedSubarray_test();

        //findPaths_test();

        //distributeCandies_test();

        //isSubtree_test();

        //checkInclusion_test();

        //matrixReshape_test();

        //arrayNesting_test();

        //nearestPalindromic_test();

        //findTilt_test();

        //arrayPairSum_test();

        //subarraySum_test();

        //ntreedepth_test();

        //QTree_test();

        //nextGreaterElement3_test();

        //leastBricks_test();

        //findMaxLength_test();

        //checkSubarraySum_test();

        // findLUSlength2_test();

        // findLUSlength_test();

        // change_test();

        // findMinMoves_test();

        // longestPalindromeSubseq_test();

        // findBottomLeftValue_test();

        // findFrequentTreeSum_test();

        // findRelativeRanks_test();

        // convertToBase7_test();

        // nextGreaterElements_test();

        // findMaximizedCapital_test();

        // findMode_test();

        // findWords_test();

        // findDiagonalOrder_test();

        // nextGreaterElement_test();

        // findMinStep_test();

        // predictTheWinner_test();

        // findMaxConsecutiveOnes_test();

        // medianSlidingWindow_test();

        // largestPalindrome_test();

        // randomlyGeneratePointsWithinACircle_test();

        // findRadius_test();

        // findMaxForm_test();

        // makesquare_test();

        // findAllConcatenatedWordsInADict_test();
        // findSubstringInWraproundString_test();
        // canIWin_test();
        // islandPerimeter_test();
        // minMoves2_test();
        // hammingDistance_test();
        // find132pattern_test();
        // findContentChildren_test();
        // fourSumCount_test();
        // minMoves_test();
        // frequencySort_test();
        // deleteNode_test();
        // SerializingAndDeserializingForBinaryTrees_test();
        // findDisappearedNumbers_test();
        // numberOfBoomerangs_test();
        // compress_test();
        // findDuplicates_test();
        // arrangeCoins_test();
        // findKthNumber_test();
        // pathSum_test();
        // findRightInterval_test();
        // eraseOverlapIntervals_test();
        // levelOrder_test();
        // countSegments_test();
        // minMutation_test();
        // countBattleships_test();
        // canPartition_test();
        // addStrings_test();
        // thirdMax_test();
        // rob_test();
        // countBits_test();
        // reverseVowels_test();
        // topKFrequent_test();
        // intersection_test();
        // SummaryRanges_test();
        // maxSumSubMatrix_test();
        // isPerfectSquare_test();
        // largestDivisibleSubset_test();
        // getSum_test();
        // kSmallestPairs_test();
        // canConstruct_test();
        // NestedInteger_test();
        // lexicalOrder_test();
        // firstUniqChar_test();
        // lastRemaining_test();
        // isRectangleCover_test();
        // decodeString_test();
        // longestSubstring_test();
        // maxRotateFunction_test();
        // findNthDigit_test();
        // removeKdigits_test();
        // canCross_test();
        // fizzBuzz_test();
        // numberOfArithmeticSlices_test();
    }
    return 0;
}
