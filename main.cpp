#include "include/createTree.h"
#include <iostream>
#include <unordered_map>

using namespace std;

struct ListNode {
    int val;
    ListNode *next;

    ListNode() : val(0), next(nullptr) {}

    ListNode(int x) : val(x), next(nullptr) {}

    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

template<typename T>
ListNode *create_nodelist(vector<T> &nums) {
    ListNode *head = nullptr;
    ListNode *curr = nullptr;
    if (nums.size() == 0)
        return nullptr;
    else
        head = new ListNode(nums[0]);
    curr = head;
    for (int i = 1; i < nums.size(); ++i) {
        curr->next = new ListNode(nums[i]);
        curr = curr->next;
    }
    return head;
}


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

namespace cutOffTree {
    int dirs[4][2] = {{-1, 0},
                      {1,  0},
                      {0,  -1},
                      {0,  1}};

    int bfs(vector<vector<int>> &forest, int sx, int sy, int tx, int ty) {
        if (sx == tx && sy == ty) {
            return 0;
        }

        int row = forest.size();
        int col = forest[0].size();
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
        vector<vector<bool>> visited(row, vector<bool>(col, false));
        pq.emplace(0, sx * col + sy);
        visited[sx][sy] = true;
        while (!pq.empty()) {
            auto[dist, loc] = pq.top();
            pq.pop();
            for (int j = 0; j < 4; ++j) {
                int nx = loc / col + dirs[j][0];
                int ny = loc % col + dirs[j][1];
                if (nx >= 0 && nx < row && ny >= 0 && ny < col) {
                    if (!visited[nx][ny] && forest[nx][ny] > 0) {
                        if (nx == tx && ny == ty) {
                            return dist + 1;
                        }
                        pq.emplace(dist + 1, nx * col + ny);
                        visited[nx][ny] = true;
                    }
                }
            }
        }
        return -1;
    }

    int cutOffTree(vector<vector<int>> &forest) {
        vector<pair<int, int>> trees;
        int row = forest.size();
        int col = forest[0].size();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (forest[i][j] > 1) {
                    trees.emplace_back(i, j);
                }
            }
        }
        sort(trees.begin(), trees.end(), [&](const pair<int, int> &a, const pair<int, int> &b) {
            return forest[a.first][a.second] < forest[b.first][b.second];
        });

        int cx = 0;
        int cy = 0;
        int ans = 0;
        for (auto &tree : trees) {
            int steps = bfs(forest, cx, cy, tree.first, tree.second);
            if (steps == -1) {
                return -1;
            }
            ans += steps;
            cx = tree.first;
            cy = tree.second;
        }
        return ans;
    }
}

void cutOffTree_test() {
    vector<vector<int>> forest;
    forest = {{1, 2, 3},
              {0, 0, 4},
              {7, 6, 5}};
    cout << cutOffTree::cutOffTree(forest) << endl;
    forest = {{1, 2, 3},
              {0, 0, 0},
              {7, 6, 5}};
    cout << cutOffTree::cutOffTree(forest) << endl;

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
    vector<int> cards;
    cards = {4, 1, 8, 7};
    cout << "true," << judgePoint24::judgePoint24(cards) << endl;
    cards = {1, 2, 1, 2};
    cout << "false," << judgePoint24::judgePoint24(cards) << endl;
}

namespace validPalindrome {
    bool is_delete;

    bool dfs(int left, int righ, string &s) {
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
                return dfs(left + 1, righ, s) || dfs(left, righ - 1, s);
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

void validPalindrome_test() {
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

namespace calPoints {
    int calPoints(vector<string> &operations) {
        vector<int> nums;
        for (auto s : operations) {
            if (s != "C" && s != "D" && s != "+") {
                nums.push_back(atoi(s.c_str()));
            } else {
                if (s == "C") {
                    nums.pop_back();
                } else if (s == "D") {
                    nums.push_back(nums.back() * 2);
                } else {
                    nums.push_back(nums[nums.size() - 1] + nums[nums.size() - 2]);
                }
            }
        }
        return accumulate(nums.begin(), nums.end(), 0);
    }
}

void calPoints_test() {
    vector<string> operations;
    operations = {"5", "2", "C", "D", "+"};
    cout << calPoints::calPoints(operations) << endl;
    operations = {"5", "-2", "4", "C", "D", "9", "+", "+"};
    cout << calPoints::calPoints(operations) << endl;
    operations = {"1"};
    cout << calPoints::calPoints(operations) << endl;
}

namespace findRedundantConnection {
    int Find(vector<int> &parent, int index) {
        if (parent[index] != index) {
            parent[index] = Find(parent, parent[index]);
        }
        return parent[index];
    }

    void Union(vector<int> &parent, int index1, int index2) {
        parent[Find(parent, index1)] = Find(parent, index2);
    }

    vector<int> findRedundantConnection(vector<vector<int>> &edges) {
        int n = edges.size();
        vector<int> parent(n + 1);
        for (int i = 1; i <= n; ++i) {
            parent[i] = i;
        }
        for (auto &edge: edges) {
            int node1 = edge[0], node2 = edge[1];
            if (Find(parent, node1) != Find(parent, node2)) {
                Union(parent, node1, node2);
            } else {
                return edge;
            }
        }
        return vector<int>{};
    }
}

void findRedundantConnection_test() {
    vector<vector<int>> edges;
    vector<int> ans;
    edges = {{1, 2},
             {1, 3},
             {2, 3}};
    ans = findRedundantConnection::findRedundantConnection(edges);
    print_vector(ans);
    edges = {{1, 2},
             {2, 3},
             {3, 4},
             {1, 4},
             {1, 5}};
    ans = findRedundantConnection::findRedundantConnection(edges);
    print_vector(ans);
}

namespace findRedundantDirectedConnection {
    struct UnionFind {
        vector<int> ancestor;

        UnionFind(int n) {
            ancestor.resize(n);
            for (int i = 0; i < n; ++i) {
                ancestor[i] = i;
            }
        }

        int find(int index) {
            return index == ancestor[index] ? index : ancestor[index] = find(ancestor[index]);
        }

        void merge(int u, int v) {
            ancestor[find(u)] = find(v);
        }
    };

    vector<int> findRedundantDirectedConnection(vector<vector<int>> &edges) {
        int n = edges.size();
        UnionFind uf = UnionFind(n + 1);
        auto parent = vector<int>(n + 1);
        for (int i = 1; i <= n; ++i) {
            parent[i] = i;
        }
        int conflict = -1;
        int cycle = -1;
        for (int i = 0; i < n; ++i) {
            auto edge = edges[i];
            int node1 = edge[0], node2 = edge[1];
            if (parent[node2] != node2) {
                conflict = i;
            } else {
                parent[node2] = node1;
                if (uf.find(node1) == uf.find(node2)) {
                    cycle = i;
                } else {
                    uf.merge(node1, node2);
                }
            }
        }
        if (conflict < 0) {
            auto redundant = vector<int>{edges[cycle][0], edges[cycle][1]};
            return redundant;
        } else {
            auto conflictEdge = edges[conflict];
            if (cycle >= 0) {
                auto redundant = vector<int>{parent[conflictEdge[1]], conflictEdge[1]};
                return redundant;
            } else {
                auto redundant = vector<int>{conflictEdge[0], conflictEdge[1]};
                return redundant;
            }
        }
    }

}

void findRedundantDirectedConnection_test() {
    vector<vector<int>> edges;
    vector<int> ans;
    edges = {{1, 2},
             {1, 3},
             {2, 3}};
    ans = findRedundantDirectedConnection::findRedundantDirectedConnection(edges);
    print_vector(ans);
    cout << "++++++++++++++" << endl;
    edges = {{1, 2},
             {2, 3},
             {3, 4},
             {4, 1},
             {1, 5}};
    ans = findRedundantDirectedConnection::findRedundantDirectedConnection(edges);
    print_vector(ans);
}

namespace repeatedStringMatch {
    int strStr(string haystack, string needle) {
        int n = haystack.size(), m = needle.size();
        if (m == 0) {
            return 0;
        }

        long long k1 = 1e9 + 7;
        long long k2 = 1337;
        srand((unsigned) time(NULL));
        long long kMod1 = rand() % k1 + k1;
        long long kMod2 = rand() % k2 + k2;

        long long hash_needle = 0;
        for (auto c : needle) {
            hash_needle = (hash_needle * kMod2 + c) % kMod1;
        }
        long long hash_haystack = 0, extra = 1;
        for (int i = 0; i < m - 1; i++) {
            hash_haystack = (hash_haystack * kMod2 + haystack[i % n]) % kMod1;
            extra = (extra * kMod2) % kMod1;
        }
        for (int i = m - 1; (i - m + 1) < n; i++) {
            hash_haystack = (hash_haystack * kMod2 + haystack[i % n]) % kMod1;
            if (hash_haystack == hash_needle) {
                return i - m + 1;
            }
            hash_haystack = (hash_haystack - extra * haystack[(i - m + 1) % n]) % kMod1;
            hash_haystack = (hash_haystack + kMod1) % kMod1;
        }
        return -1;
    }

    int repeatedStringMatch(string a, string b) {
        int an = a.size(), bn = b.size();
        int index = strStr(a, b);
        if (index == -1) {
            return -1;
        }
        if (an - index >= bn) {
            return 1;
        }
        return (bn + index - an - 1) / an + 2;
    }
}

void repeatedStringMatch_test() {
    string a, b;
    a = "abcd";
    b = "cdabcdab";
    cout << repeatedStringMatch::repeatedStringMatch(a, b) << endl;
    a = "a";
    b = "aa";
    cout << repeatedStringMatch::repeatedStringMatch(a, b) << endl;
    a = "a";
    b = "a";
    cout << repeatedStringMatch::repeatedStringMatch(a, b) << endl;
    a = "abc";
    b = "wxyz";
    cout << repeatedStringMatch::repeatedStringMatch(a, b) << endl;
}

namespace maxSumOfThreeSubarrays {
    vector<int> maxSumOfThreeSubarrays(vector<int> &nums, int k) {
        vector<int> ans;
        int sum1 = 0, maxSum1 = 0, maxSum1Idx = 0;
        int sum2 = 0, maxSum12 = 0, maxSum12Idx1 = 0, maxSum12Idx2 = 0;
        int sum3 = 0, maxTotal = 0;
        for (int i = k * 2; i < nums.size(); ++i) {
            sum1 += nums[i - k * 2];
            sum2 += nums[i - k];
            sum3 += nums[i];
            if (i >= k * 3 - 1) {
                if (sum1 > maxSum1) {
                    maxSum1 = sum1;
                    maxSum1Idx = i - k * 3 + 1;
                }
                if (maxSum1 + sum2 > maxSum12) {
                    maxSum12 = maxSum1 + sum2;
                    maxSum12Idx1 = maxSum1Idx;
                    maxSum12Idx2 = i - k * 2 + 1;
                }
                if (maxSum12 + sum3 > maxTotal) {
                    maxTotal = maxSum12 + sum3;
                    ans = {maxSum12Idx1, maxSum12Idx2, i - k + 1};
                }
                sum1 -= nums[i - k * 3 + 1];
                sum2 -= nums[i - k * 2 + 1];
                sum3 -= nums[i - k + 1];
            }
        }
        return ans;
    }
}

void maxSumOfThreeSubarrays_test() {
    vector<int> nums, ans;
    int k;
    nums = {1, 2, 1, 2, 6, 7, 5, 1};
    k = 2;
    ans = maxSumOfThreeSubarrays::maxSumOfThreeSubarrays(nums, k);
    print_vector(ans);
    nums = {1, 2, 1, 2, 1, 2, 1, 2, 1};
    k = 2;
    ans = maxSumOfThreeSubarrays::maxSumOfThreeSubarrays(nums, k);
    print_vector(ans);
}

namespace minStickers {
    int minStickers(vector<string> &stickers, string target) {
        int m = target.size();//首先获取target的长度，记m,mask为二进制表示的target的每一位
        //给定一个1<<m大小的vector每个位置来表示当选用mask为当前值时候（这个可能没表述好，看官方解答
        //所对应的最少的stricker的个数
        //初值给-1，方便下面判断
        vector<int> dp(1 << m, -1);
        dp[0] = 0;//表示空的子串需要0个sticker就可以完成
        //function c++11,表示一个函数指针？（我也刚刚看不确定)表示类模板，模板参数为int(int),外层int表示返回值，内层
        //int表示参数类型，用lambda来构建该函数
        function<int(int)> helper = [&](int mask) {
            if (dp[mask] != -1) {
                return dp[mask];//记忆化搜索需要，如果已经算到过当前mask所对应的所需要的最小stricker数，直接返回
                //当然，以案例 strickers="a" ,target="aa"来说，需要"a"2次，当搜索到最后一轮的时候，初值起作用，直接返回0
            }
            //关于这个，我尝试了一下给一个比较大的数都可以，因为最多拼接次数可以为当前字符串的长度m
            //m+1之后dp在每次min求解时会被更新，min求解见43行，
            dp[mask] = m + 1;
            for (auto &stricker:stickers) {//既然要求出最少的次数，那么对每个stricker都需要进行求解
                //mask为传入这个函数的时候所拥有的，对于target来说，也就是当前还剩下的东西，因为对于每个stricker求解都需要
                //mask作为初始状态，为了不改变这个初始状态，再给一个参数left，这个left也会用在之后的递归中，
                //不给那就没办法作为递归参数传给下一层的mask了
                int left = mask;
                vector<int> count(26, 0);//没啥好说的，（总不可能是24个字母吧 :)
                for (auto &cr:stricker) {
                    count[cr - 'a']++;//记录当前的stricker有多少字母组成，每个字母的个数
                }
                //有了上述条件之后对left可以进行求解了,m已经在开头说了，表示为target长度，那么对整个长度里所有的字母进行遍历呗
                for (int i = 0; i < m; i++) {
                    //条件是mask里当前位置的字母存在，并且stricker里对应的字母出现次数还够
                    if (((mask >> i) & 1) && count[target[i] - 'a'] > 0) {
                        count[target[i] - 'a']--;//用掉了嘛，那就--
                        left ^= 1 << i;//left在之前和mask一样，只是left用来表示被裁剪之后，因此某个字母被剪掉了
                        //可以用异或运算，相同为0，相异为1，当前字母在mask里存在，为1，1<<i表示该字母的位置，异或一下，就为0了
                        //这时候left可以表示除掉该字母之后的情况
                    }
                    //经过上述的操作，left和mask的区别是left已经除去了16行开始的循环里的某个stricker里的字母了，但是
                    //可能还没全部归0，因此要进行递归求解，在求解之前先判断当前的剩下的left是否小于mask，没变化的话就不需要继续递归啦，
                    //因为stricker里面完全没有对应的可以删去的字符，那就直接下一个
                }
                if (left < mask) {
                    //求个最小，要么当前状态，要么剩下的字母去求解，得到的stricker的次数再加上当前次，也就是+1
                    dp[mask] = min(dp[mask], helper(left) + 1);
                }
            }
            return dp[mask];//返回目标需要次数
        };//lambda函数是个表达式，得加;
        //一开始传入多少呢？以target="aa"来说，m=2, 1<<2=4,(1<<2)-1=3,mask=3,换成二进制也就是11表示的是当前aa这两位的情况
        //同理对于其他的，比如target="aaaaa" m=5, 1<<5=32,(1<<2)-1=31,二进制表示11111，表示所有都存在的情况
        int ans = helper((1 << m) - 1);
        return ans > m ? -1 : ans;
    }
}

void minStickers_test() {
    vector<string> stickers;
    string target;
    stickers = {"with", "example", "science"};
    target = "thehat";
    cout << minStickers::minStickers(stickers, target) << endl;
    stickers = {"notice", "possible"};
    target = "basicbasic";
    cout << minStickers::minStickers(stickers, target) << endl;
}

namespace topKFrequent {
    vector<string> topKFrequent(vector<string> &words, int k) {
        unordered_map<string, int> cnt;
        for (auto &word : words) {
            ++cnt[word];
        }
        vector<string> rec;
        for (auto&[key, value] : cnt) {
            rec.emplace_back(key);
        }
        sort(rec.begin(), rec.end(), [&](const string &a, const string &b) -> bool {
            return cnt[a] == cnt[b] ? a < b : cnt[a] > cnt[b];
        });
        rec.erase(rec.begin() + k, rec.end());
        return rec;
    }
}

void topKFrequent692_test() {
    vector<string> words, ans;
    int k;
    words = {"i", "love", "leetcode", "i", "love", "coding"};
    k = 2;
    ans = topKFrequent::topKFrequent(words, k);
    print_vector(ans);
    words = {"the", "day", "is", "sunny", "the", "the", "the", "sunny", "is", "is"};
    k = 4;
    ans = topKFrequent::topKFrequent(words, k);
    print_vector(ans);
}

namespace hasAlternatingBits {
    bool hasAlternatingBits(int n) {
        long a = n ^(n >> 1);
        return (a & (a + 1)) == 0;
    }
}

void hasAlternatingBits_test() {
    int n;
    n = 5;
    cout << hasAlternatingBits::hasAlternatingBits(n) << endl;
    n = 7;
    cout << hasAlternatingBits::hasAlternatingBits(n) << endl;
    n = 11;
    cout << hasAlternatingBits::hasAlternatingBits(n) << endl;
}

namespace maxAreaOfIslan {
    int dfs(vector<vector<int>> &grid, int cur_i, int cur_j) {
        if (cur_i < 0 || cur_j < 0 || cur_i == grid.size() || cur_j == grid[0].size() || grid[cur_i][cur_j] != 1) {
            return 0;
        }
        grid[cur_i][cur_j] = 0;
        int di[4] = {0, 0, 1, -1};
        int dj[4] = {1, -1, 0, 0};
        int ans = 1;
        for (int index = 0; index != 4; ++index) {
            int next_i = cur_i + di[index], next_j = cur_j + dj[index];
            ans += dfs(grid, next_i, next_j);
        }
        return ans;
    }

    int maxAreaOfIsland(vector<vector<int>> &grid) {
        int ans = 0;
        for (int i = 0; i != grid.size(); ++i) {
            for (int j = 0; j != grid[0].size(); ++j) {
                ans = max(ans, dfs(grid, i, j));
            }
        }
        return ans;
    }
}

void maxAreaOfIsland_test() {
    vector<vector<int>> grid;
    grid = {{0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0},
            {0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0},
            {0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0}};
    cout << maxAreaOfIslan::maxAreaOfIsland(grid) << endl;
    grid = {{0, 0, 0, 0, 0, 0, 0, 0}};
    cout << maxAreaOfIslan::maxAreaOfIsland(grid) << endl;
}

namespace countBinarySubstrings {
    int countBinarySubstrings(string s) {
        int ptr = 0, n = s.size(), last = 0, ans = 0;
        while (ptr < n) {
            char c = s[ptr];
            int count = 0;
            while (ptr < n && s[ptr] == c) {
                ++ptr;
                ++count;
            }
            ans += min(count, last);
            last = count;
        }
        return ans;
    }
}

void countBinarySubstrings_test() {
    string s;
    s = "00110011";
    cout << countBinarySubstrings::countBinarySubstrings(s) << endl;
    s = "10101";
    cout << countBinarySubstrings::countBinarySubstrings(s) << endl;
}

namespace findShortestSubArray {
    int findShortestSubArray(vector<int> &nums) {
        unordered_map<int, vector<int>> mp;
        int n = nums.size();
        for (int i = 0; i < n; i++) {
            if (mp.count(nums[i])) {
                mp[nums[i]][0]++;
                mp[nums[i]][2] = i;
            } else {
                mp[nums[i]] = {1, i, i};
            }
        }
        int maxNum = 0, minLen = 0;
        for (auto&[_, vec] : mp) {
            if (maxNum < vec[0]) {
                maxNum = vec[0];
                minLen = vec[2] - vec[1] + 1;
            } else if (maxNum == vec[0]) {
                if (minLen > vec[2] - vec[1] + 1) {
                    minLen = vec[2] - vec[1] + 1;
                }
            }
        }
        return minLen;
    }
}

void findShortestSubArray_test() {
    vector<int> nums;
    nums = {1, 2, 2, 3, 1};
    cout << findShortestSubArray::findShortestSubArray(nums) << endl;
    nums = {1, 2, 2, 3, 1, 4, 2};
    cout << findShortestSubArray::findShortestSubArray(nums) << endl;
}

namespace canPartitionKSubsets {
    bool canPartitionKSubsets(vector<int> &nums, int k) {
        int all = accumulate(nums.begin(), nums.end(), 0);
        if (all % k > 0) {
            return false;
        }
        int per = all / k;
        sort(nums.begin(), nums.end());
        if (nums.back() > per) {
            return false;
        }
        int n = nums.size();
        vector<bool> dp(1 << n, false);
        vector<int> curSum(1 << n, 0);
        dp[0] = true;
        auto ss = 1 << n;
        cout << ss << endl;
        for (int i = 0; i < 1 << n; i++) {
            if (!dp[i]) {
                continue;
            }
            for (int j = 0; j < n; j++) {
                if (curSum[i] + nums[j] > per) {
                    break;
                }
                if (((i >> j) & 1) == 0) {
                    int next = i | (1 << j);
                    if (!dp[next]) {
                        curSum[next] = (curSum[i] + nums[j]) % per;
                        dp[next] = true;
                    }
                }
            }
        }
        return dp[(1 << n) - 1];
    }
}

void canPartitionKSubsets_test() {
    vector<int> nums;
    int k;
    nums = {4, 3, 2, 3, 5, 2, 1};
    k = 4;
    cout << canPartitionKSubsets::canPartitionKSubsets(nums, k) << endl;
    nums = {1, 2, 3, 4};
    k = 3;
    cout << canPartitionKSubsets::canPartitionKSubsets(nums, k) << endl;
}

namespace fallingSquares {
    vector<int> fallingSquares(vector<vector<int>> &positions) {
        int n = positions.size();
        vector<int> ret(n);
        map<int, int> heightMap;
        heightMap[0] = 0; // 初始时从 0 开始的所有点的堆叠高度都是 0
        for (int i = 0; i < n; i++) {
            int size = positions[i][1];
            int left = positions[i][0], right = positions[i][0] + positions[i][1] - 1;
            auto lp = heightMap.upper_bound(left), rp = heightMap.upper_bound(right);

            //int rHeight = prev(rp)->second; // 记录 right + 1 对应的堆叠高度（如果 right + 1 不在 heightMap 中）
            int rHeight = rp->second;

            // 更新第 i 个掉落的方块的堆叠高度
            int height = 0;
            for (auto p = prev(lp); p != rp; p++) {
                height = max(height, p->second + size);
            }

            // 清除 heightMap 中位于 (left, right] 内的点
            heightMap.erase(lp, rp);

            heightMap[left] = height; // 更新 left 的变化
            if (rp == heightMap.end() || rp->first != right + 1) { // 如果 right + 1 不在 heightMap 中，更新 right + 1 的变化
                heightMap[right + 1] = rHeight;
            }
            ret[i] = i > 0 ? max(ret[i - 1], height) : height;
        }
        return ret;
    }
}

void fallingSquares_test() {
    vector<vector<int>> positions;
    vector<int> ans;
    positions = {{1, 2},
                 {2, 3},
                 {6, 1}};
    ans = fallingSquares::fallingSquares(positions);
    print_vector(ans);
    positions = {{100, 100},
                 {200, 100}};
    ans = fallingSquares::fallingSquares(positions);
    print_vector(ans);

}

namespace searchBST {
    TreeNode::TreeNode *searchBST(TreeNode::TreeNode *root, int val) {
        if (root == nullptr)
            return nullptr;
        if (root->val == val)
            return root;
        return searchBST(val < root->val ? root->left : root->right, val);
    }
}

void searchBST_test() {
    TreeNode::TreeNode *ans, *root;
    vector<int> root_val;
    int val;
    root_val = {4, 2, 7, 1, 3};
    val = 2;
    root = create_treenode(root_val);
    ans = searchBST::searchBST(root, val);
    cout << TreeNode::print_tree(ans) << endl;
    root_val = {4, 2, 7, 1, 3};
    val = 5;
    root = create_treenode(root_val);
    ans = searchBST::searchBST(root, val);
    cout << TreeNode::print_tree(ans) << endl;
}

namespace insertIntoBST {
    TreeNode::TreeNode *insertIntoBST(TreeNode::TreeNode *root, int val) {
        if (root == nullptr) {
            return new TreeNode::TreeNode(val);
        }
        TreeNode::TreeNode *pos = root;
        while (pos != nullptr) {
            if (val < pos->val) {
                if (pos->left == nullptr) {
                    pos->left = new TreeNode::TreeNode(val);
                    break;
                } else {
                    pos = pos->left;
                }
            } else {
                if (pos->right == nullptr) {
                    pos->right = new TreeNode::TreeNode(val);
                    break;
                } else {
                    pos = pos->right;
                }
            }
        }
        return root;
    }
}

void insertIntoBST_test() {
    vector<int> vals;
    int val;
    TreeNode::TreeNode *root, *ans;
    vals = {4, 2, 7, 1, 3};
    val = 5;
    root = create_treenode(vals);
    ans = insertIntoBST::insertIntoBST(root, val);
    cout << TreeNode::print_tree(ans) << endl;
    vals = {40, 20, 60, 10, 30, 50, 70};
    val = 25;
    root = create_treenode(vals);
    ans = insertIntoBST::insertIntoBST(root, val);
    cout << TreeNode::print_tree(ans) << endl;
    vals = {4, 2, 7, 1, 3, 0, 0, 0, 0, 0, 0};
    val = 25;
    root = create_treenode(vals);
    ans = insertIntoBST::insertIntoBST(root, val);
    cout << TreeNode::print_tree(ans) << endl;
}

namespace binarySearch {
    int binarySearch(vector<int> &nums, int target, int left, int righ) {
        if (left > righ || (left == righ && nums[left] != target)) {
            return -1;
        }
        int mid = left + (righ - left) / 2;
        if (nums[mid] > target) {
            return binarySearch(nums, target, left, mid);
        } else if (nums[mid] < target) {
            return binarySearch(nums, target, mid + 1, righ);
        } else {
            return mid;
        }
    }

    int search(vector<int> &nums, int target) {
        return binarySearch(nums, target, 0, nums.size() - 1);
    }
}

void binarySearch_test() {
    vector<int> nums;
    int target;
    nums = {-1, 0, 3, 5, 9, 12};
    target = 9;
    cout << binarySearch::search(nums, target) << endl;
    nums = {-1, 0, 3, 5, 9, 12};
    target = 2;
    cout << binarySearch::search(nums, target) << endl;
}

namespace minimumDeleteSum {
    int minimumDeleteSum(string s1, string s2) {
        int m = s1.size();
        int n = s2.size();
        vector<vector<int>> dp(m + 1, vector<int>(n + 1));

        for (int i = 1; i <= m; ++i) {
            dp[i][0] = dp[i - 1][0] + s1[i - 1];
        }
        for (int j = 1; j <= n; ++j) {
            dp[0][j] = dp[0][j - 1] + s2[j - 1];
        }
        for (int i = 1; i <= m; i++) {
            char c1 = s1[i - 1];
            for (int j = 1; j <= n; j++) {
                char c2 = s2[j - 1];
                if (c1 == c2) {
                    dp[i][j] = dp[i - 1][j - 1];
                } else {
                    dp[i][j] = min(dp[i - 1][j] + s1[i - 1], dp[i][j - 1] + s2[j - 1]);
                }
            }
        }

        return dp[m][n];
    }
}

void minimumDeleteSum_test() {
    string s1, s2;
    s1 = "sea";
    s2 = "eat";
    cout << minimumDeleteSum::minimumDeleteSum(s1, s2) << endl;
    s1 = "delete";
    s2 = "leet";
    cout << minimumDeleteSum::minimumDeleteSum(s1, s2) << endl;
}

namespace numSubarrayProductLessThanK {
    int numSubarrayProductLessThanK(vector<int> &nums, int k) {
        int n = nums.size(), ans = 0;
        int prod = 1, i = 0;
        for (int j = 0; j < n; j++) {
            prod *= nums[j];
            while (i <= j && prod >= k) {
                prod /= nums[i];
                i++;
            }
            ans += j - i + 1;
        }
        return ans;
    }
}

void numSubarrayProductLessThanK_test() {
    vector<int> nums;
    int k;
    nums = {10, 5, 2, 6};
    k = 100;
    cout << numSubarrayProductLessThanK::numSubarrayProductLessThanK(nums, k) << endl;
    nums = {1, 2, 3};
    k = 0;
    cout << numSubarrayProductLessThanK::numSubarrayProductLessThanK(nums, k) << endl;
}

namespace maxProfit {
    int maxProfit(vector<int> &prices, int fee) {
        int n = prices.size();
        vector<vector<int>> dp(n, vector<int>(2));
        dp[0][0] = 0, dp[0][1] = -prices[0];
        for (int i = 1; i < n; ++i) {
            dp[i][0] = max(dp[i - 1][0], dp[i - 1][1] + prices[i] - fee);
            dp[i][1] = max(dp[i - 1][1], dp[i - 1][0] - prices[i]);
        }
        return dp[n - 1][0];
    }
}

void maxProfit_test() {
    vector<int> prices;
    int fee;
    prices = {1, 3, 2, 8, 4, 9};
    fee = 2;
    cout << maxProfit::maxProfit(prices, fee) << endl;
    prices = {1, 3, 7, 5, 10, 3};
    fee = 3;
    cout << maxProfit::maxProfit(prices, fee) << endl;
}

namespace findLength {
    int findLength(vector<int> &nums1, vector<int> &nums2) {
        int n = nums1.size();
        int m = nums2.size();
        int ans = 0;
        vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
        for (int i = n - 1; i >= 0; --i) {
            for (int j = m - 1; j >= 0; --j) {
                dp[i][j] = nums1[i] == nums2[j] ? dp[i + 1][j + 1] + 1 : 0;
                ans = max(ans, dp[i][j]);
            }
        }
        return ans;
    }
}

void findLength_test() {
    vector<int> nums1, nums2;
    nums1 = {1, 2, 3, 2, 1};
    nums2 = {3, 2, 1, 4, 7};
    cout << findLength::findLength(nums1, nums2) << endl;
    nums1 = {0, 0, 0, 0, 0};
    nums2 = {0, 0, 0, 0, 0};
    cout << findLength::findLength(nums1, nums2) << endl;
}

namespace smallestDistancePair {
    int smallestDistancePair(vector<int> &nums, int k) {
        sort(nums.begin(), nums.end());
        int n = nums.size(), left = 0, right = nums.back() - nums.front();
        while (left <= right) {
            int mid = (left + right) / 2;
            int cnt = 0;
            for (int i = 0, j = 0; j < n; j++) {
                while (nums[j] - nums[i] > mid) {
                    i++;
                }
                cnt += j - i;
            }
            if (cnt >= k) {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        return left;
    }
}

void smallestDistancePair_test() {
    vector<int> nums;
    int k;
    nums = {62, 100, 4};
    k = 2;
    cout << smallestDistancePair::smallestDistancePair(nums, k) << endl;
}

namespace longestWord {
    class Trie {
    public:
        Trie() {
            this->children = vector<Trie *>(26, nullptr);
            this->is_end = false;
        }

        bool insert(const string &word) {
            Trie *node = this;
            for (const auto &ch : word) {
                int index = ch - 'a';
                if (node->children[index] == nullptr) {
                    node->children[index] = new Trie();
                }
                node = node->children[index];
            }
            node->is_end = true;
            return true;
        }

        bool search(const string &word) {
            Trie *node = this;
            for (const auto &ch : word) {
                int index = ch - 'a';
                if (node->children[index] == nullptr || !node->children[index]->is_end) {
                    return false;
                }
                node = node->children[index];
            }
            return node != nullptr && node->is_end;
        }

    private:
        vector<Trie *> children;
        bool is_end;
    };

    string longestWord(vector<string> &words) {
        Trie trie;
        for (const auto &word : words) {
            trie.insert(word);
        }
        string longest = "";
        for (const auto &word : words) {
            if (trie.search(word)) {
                if (word.size() > longest.size() || (word.size() == longest.size() && word < longest)) {
                    longest = word;
                }
            }
        }
        return longest;
    }
}

void longestWord_test() {
    vector<string> words;
    string ans;
    words = {"w", "wo", "wor", "worl", "world"};
    ans = longestWord::longestWord(words);
    cout << ans << endl;
    words = {"a", "banana", "app", "appl", "ap", "apply", "apple"};
    ans = longestWord::longestWord(words);
    cout << ans << endl;
}

namespace accountsMerge {
    class UnionFind {
    public:
        vector<int> parent;

        UnionFind(int n) {
            parent.resize(n);
            for (int i = 0; i < n; i++) {
                parent[i] = i;
            }
        }

        void unionSet(int index1, int index2) {
            parent[find(index2)] = find(index1);
        }

        int find(int index) {
            if (parent[index] != index) {
                parent[index] = find(parent[index]);
            }
            return parent[index];
        }
    };

    vector<vector<string>> accountsMerge(vector<vector<string>> &accounts) {
        map<string, int> emailToIndex;
        map<string, string> emailToName;
        int emailsCount = 0;
        for (auto &account : accounts) {
            string &name = account[0];
            int size = account.size();
            for (int i = 1; i < size; i++) {
                string &email = account[i];
                if (!emailToIndex.count(email)) {
                    emailToIndex[email] = emailsCount++;
                    emailToName[email] = name;
                }
            }
        }
        UnionFind uf(emailsCount);
        for (auto &account : accounts) {
            string &firstEmail = account[1];
            int firstIndex = emailToIndex[firstEmail];
            int size = account.size();
            for (int i = 2; i < size; i++) {
                string &nextEmail = account[i];
                int nextIndex = emailToIndex[nextEmail];
                uf.unionSet(firstIndex, nextIndex);
            }
        }
        map<int, vector<string>> indexToEmails;
        for (auto&[email, _] : emailToIndex) {
            int index = uf.find(emailToIndex[email]);
            vector<string> &account = indexToEmails[index];
            account.emplace_back(email);
            indexToEmails[index] = account;
        }
        vector<vector<string>> merged;
        for (auto&[_, emails] : indexToEmails) {
            sort(emails.begin(), emails.end());
            string &name = emailToName[emails[0]];
            vector<string> account;
            account.emplace_back(name);
            for (auto &email : emails) {
                account.emplace_back(email);
            }
            merged.emplace_back(account);
        }
        return merged;
    }
}

void accountsMerge_test() {
    vector<vector<string>> accounts, ans;
    accounts = {{"John", "johnsmith@mail.com", "john00@mail.com"},
                {"John", "johnnybravo@mail.com"},
                {"John", "johnsmith@mail.com", "john_newyork@mail.com"},
                {"Mary", "mary@mail.com"}};
    ans = accountsMerge::accountsMerge(accounts);
    for (auto &word : ans) {
        print_vector(word);
    }
    cout << "--------------" << endl;
    accounts = {{"John", "john00@mail.com", "john_newyork@mail.com", "johnsmith@mail.com"},
                {"John", "johnnybravo@mail.com"},
                {"Mary", "mary@mail.com"}};
    ans = accountsMerge::accountsMerge(accounts);
    for (auto &word : ans) {
        print_vector(word);
    }
    cout << "--------------" << endl;
}

namespace pivotIndex {
    int pivotIndex(vector<int> &nums) {
        if (nums.size() == 1) {
            return 0;
        }
        if (nums.size() == 0)
            return -1;
        int n = nums.size();
        int leftsum = 0, righsum = accumulate(nums.begin() + 1, nums.end(), 0.0);
        for (int i = 0; i < n - 1; ++i) {
            if (leftsum == righsum) {
                return i;
            } else {
                leftsum += nums[i];
                righsum -= nums[i + 1];
            }
        }
        if (leftsum == righsum) {
            return n - 1;
        }
        return -1;
    }
}

void pivotIndex_test() {
    vector<int> nums;
    nums = {1, 7, 3, 6, 5, 6};
    cout << pivotIndex::pivotIndex(nums) << endl;
    nums = {1, 2, 3};
    cout << pivotIndex::pivotIndex(nums) << endl;
    nums = {2, -1, 1};
    cout << pivotIndex::pivotIndex(nums) << endl;
    nums = {-1, 1, 2};
    cout << pivotIndex::pivotIndex(nums) << endl;
}

namespace splitListToParts {
    vector<ListNode *> splitListToParts(ListNode *head, int k) {
        int n = 0;
        ListNode *temp = head;
        while (temp != nullptr) {
            n++;
            temp = temp->next;
        }
        int quotient = n / k, remainder = n % k;

        vector<ListNode *> parts(k, nullptr);
        ListNode *curr = head;
        for (int i = 0; i < k && curr != nullptr; i++) {
            parts[i] = curr;
            int partSize = quotient + (i < remainder ? 1 : 0);
            for (int j = 1; j < partSize; j++) {
                curr = curr->next;
            }
            ListNode *next = curr->next;
            curr->next = nullptr;
            curr = next;
        }
        return parts;
    }

    void print(vector<ListNode *> nodelist) {
        for (auto list : nodelist) {
            cout << "{";
            while (list != nullptr) {
                cout << list->val;
                list = list->next;
                if (list) {
                    cout << ",";
                }
            }
            cout << "} " << endl;
        }
    }
}

void splitListToParts_test() {
    vector<int> nums;
    int k;
    vector<ListNode *> ans;
    nums = {1, 2, 3};
    k = 5;
    ListNode *head;
    head = create_nodelist(nums);
    ans = splitListToParts::splitListToParts(head, k);
    splitListToParts::print(ans);
    cout << "+++++++++++" << endl;
    nums = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    k = 3;
    head = create_nodelist(nums);
    ans = splitListToParts::splitListToParts(head, k);
    splitListToParts::print(ans);
    cout << "+++++++++++" << endl;
}

namespace countOfAtoms {
    string countOfAtoms(string formula) {
        int i = 0, n = formula.length();

        auto parseAtom = [&]() -> string {
            string atom;
            atom += formula[i++]; // 扫描首字母
            while (i < n && islower(formula[i])) {
                atom += formula[i++]; // 扫描首字母后的小写字母
            }
            return atom;
        };

        auto parseNum = [&]() -> int {
            if (i == n || !isdigit(formula[i])) {
                return 1; // 不是数字，视作 1
            }
            int num = 0;
            while (i < n && isdigit(formula[i])) {
                num = num * 10 + int(formula[i++] - '0'); // 扫描数字
            }
            return num;
        };

        stack<unordered_map<string, int>> stk;
        stk.push({});
        while (i < n) {
            char ch = formula[i];
            if (ch == '(') {
                i++;
                stk.push({}); // 将一个空的哈希表压入栈中，准备统计括号内的原子数量
            } else if (ch == ')') {
                i++;
                int num = parseNum(); // 括号右侧数字
                auto atomNum = stk.top();
                stk.pop(); // 弹出括号内的原子数量
                for (auto &[atom, v] : atomNum) {
                    stk.top()[atom] += v * num; // 将括号内的原子数量乘上 num，加到上一层的原子数量中
                }
            } else {
                string atom = parseAtom();
                int num = parseNum();
                stk.top()[atom] += num; // 统计原子数量
            }
        }

        auto &atomNum = stk.top();
        vector<pair<string, int>> pairs;
        for (auto &[atom, v] : atomNum) {
            pairs.emplace_back(atom, v);
        }
        sort(pairs.begin(), pairs.end());

        string ans;
        for (auto &p : pairs) {
            ans += p.first;
            if (p.second > 1) {
                ans += to_string(p.second);
            }
        }
        return ans;
    }
}

void countOfAtoms_test() {
    string formula;
    formula = "H2O";
    cout << countOfAtoms::countOfAtoms(formula) << endl;
    formula = "Mg(OH)2";
    cout << countOfAtoms::countOfAtoms(formula) << endl;
    formula = "K4(ON(SO3)2)2";
    cout << countOfAtoms::countOfAtoms(formula) << endl;
}

namespace selfDividingNumbers {
    bool isSelfDividing(int num) {
        int temp = num;
        while (temp > 0) {
            int digit = temp % 10;
            if (digit == 0 || num % digit != 0) {
                return false;
            }
            temp /= 10;
        }
        return true;
    }

    vector<int> selfDividingNumbers(int left, int right) {
        vector<int> ans;
        for (int i = left; i <= right; i++) {
            if (isSelfDividing(i)) {
                ans.emplace_back(i);
            }
        }
        return ans;
    }
}

void selfDividingNumbers_test() {
    int left, right;
    vector<int> ans;
    left = 1, right = 22;
    ans = selfDividingNumbers::selfDividingNumbers(left, right);
    print_vector(ans);
    left = 47, right = 85;
    ans = selfDividingNumbers::selfDividingNumbers(left, right);
    print_vector(ans);
}

namespace MyCalendar {
    class MyCalendar {
        set<pair<int, int>> booked;

    public:
        bool book(int start, int end) {
            auto it = booked.lower_bound({end, 0});
            if (it == booked.begin() || (--it)->second <= start) {
                booked.emplace(start, end);
                return true;
            }
            return false;
        }
    };
}

void MyCalendar_test() {
    MyCalendar::MyCalendar test;
    cout << test.book(10, 20) << endl;
    cout << test.book(15, 25) << endl;
    cout << test.book(20, 30) << endl;
    cout << "+++++++++" << endl;
}

namespace countPalindromicSubsequences {
    const int MOD = 1E9 + 7;

    int countPalindromicSubsequences(string s) {
        int n = s.size();
        vector<vector<vector<int>>> dp(4, vector<vector<int>>(n, vector<int>(n, 0)));
        for (int i = 0; i < n; i++) {
            dp[s[i] - 'a'][i][i] = 1;
        }

        for (int len = 2; len <= n; len++) {
            for (int i = 0, j = len - 1; j < n; i++, j++) {
                for (char c = 'a', k = 0; c <= 'd'; c++, k++) {
                    if (s[i] == c && s[j] == c) {
                        dp[k][i][j] = (2LL + dp[0][i + 1][j - 1] + dp[1][i + 1][j - 1] + dp[2][i + 1][j - 1] +
                                       dp[3][i + 1][j - 1]) % MOD;
                    } else if (s[i] == c) {
                        dp[k][i][j] = dp[k][i][j - 1];
                    } else if (s[j] == c) {
                        dp[k][i][j] = dp[k][i + 1][j];
                    } else {
                        dp[k][i][j] = dp[k][i + 1][j - 1];
                    }
                }
            }
        }

        int res = 0;
        for (int i = 0; i < 4; i++) {
            res = (res + dp[i][0][n - 1]) % MOD;
        }
        return res;
    }
}

void countPalindromicSubsequences_test() {
    string s;
    s = "bccb";
    cout << countPalindromicSubsequences::countPalindromicSubsequences(s) << endl;
    s = "abcdabcdabcdabcdabcdabcdabcdabcddcbadcbadcbadcbadcbadcbadcbadcba";
    cout << countPalindromicSubsequences::countPalindromicSubsequences(s) << endl;
}

namespace floodFill {
    const int dx[4] = {1, 0, 0, -1};
    const int dy[4] = {0, 1, -1, 0};

    void dfs(vector<vector<int>> &image, int x, int y, int currColor, int color) {
        if (image[x][y] == currColor) {
            image[x][y] = color;
            for (int i = 0; i < 4; i++) {
                int mx = x + dx[i], my = y + dy[i];
                if (mx >= 0 && mx < image.size() && my >= 0 && my < image[0].size()) {
                    dfs(image, mx, my, currColor, color);
                }
            }
        }
    }

    vector<vector<int>> floodFill(vector<vector<int>> &image, int sr, int sc, int color) {
        int currColor = image[sr][sc];
        if (currColor != color) {
            dfs(image, sr, sc, currColor, color);
        }
        return image;
    }
}

void floodFill_test() {
    vector<vector<int>> image, ans;
    int sr, sc, color;
    image = {{1, 1, 1},
             {1, 1, 0},
             {1, 0, 1}}, sr = 1, sc = 1, color = 2;
    ans = floodFill::floodFill(image, sr, sc, color);
    for (auto line : ans) {
        print_vector(line);
    }
}

namespace asteroidCollision {
    vector<int> asteroidCollision(vector<int> &asteroids) {
        vector<int> st;
        for (auto aster : asteroids) {
            bool alive = true;
            while (alive && aster < 0 && !st.empty() && st.back() > 0) {
                alive = st.back() < -aster; // aster 是否存在
                if (st.back() <= -aster) {  // 栈顶行星爆炸
                    st.pop_back();
                }
            }
            if (alive) {
                st.push_back(aster);
            }
        }
        return st;
    }
}

void asteroidCollision_test() {
    vector<int> asteroids, ans;
    asteroids = {-2, -1, 1, 2};
    ans = asteroidCollision::asteroidCollision(asteroids);
    print_vector(ans);
    asteroids = {5, 10, -5};
    ans = asteroidCollision::asteroidCollision(asteroids);
    print_vector(ans);
    asteroids = {8, -8};
    ans = asteroidCollision::asteroidCollision(asteroids);
    print_vector(ans);
}

namespace dailyTemperatures {
    vector<int> dailyTemperatures(vector<int> &temperatures) {
        int n = temperatures.size();
        vector<int> ans(n);
        stack<int> s;
        for (int i = 0; i < n; ++i) {
            while (!s.empty() && temperatures[i] > temperatures[s.top()]) {
                int previousIndex = s.top();
                ans[previousIndex] = i - previousIndex;
                s.pop();
            }
            s.push(i);
        }
        return ans;
    }
}

void dailyTemperatures_test() {
    vector<int> temperatures, ans;
    temperatures = {73, 74, 75, 71, 69, 72, 76, 73};
    ans = dailyTemperatures::dailyTemperatures(temperatures);
    print_vector(ans);
    temperatures = {30, 40, 50, 60};
    ans = dailyTemperatures::dailyTemperatures(temperatures);
    print_vector(ans);
}

namespace deleteAndEarn {
    int rob(vector<int> &nums) {
        int size = nums.size();
        if (size == 1) {
            return nums[0];
        }
        int first = nums[0], second = max(nums[0], nums[1]);
        for (int i = 2; i < size; i++) {
            int temp = second;
            second = max(first + nums[i], second);
            first = temp;
        }
        return second;
    }

    int deleteAndEarn(vector<int> &nums) {
        int n = nums.size();
        int ans = 0;
        sort(nums.begin(), nums.end());
        vector<int> sum = {nums[0]};
        for (int i = 1; i < n; ++i) {
            int val = nums[i];
            if (val == nums[i - 1]) {
                sum.back() += val;
            } else if (val == nums[i - 1] + 1) {
                sum.push_back(val);
            } else {
                ans += rob(sum);
                sum = {val};
            }
        }
        ans += rob(sum);
        return ans;
    }
}

void deleteAndEarn_test() {
    vector<int> nums;
    nums = {3, 4, 2};
    cout << deleteAndEarn::deleteAndEarn(nums) << endl;
    nums = {2, 2, 3, 3, 3, 4};
    cout << deleteAndEarn::deleteAndEarn(nums) << endl;
}

namespace cherryPickup {
    int cherryPickup(vector<vector<int>> &grid) {
        int n = grid.size();
        vector<vector<vector<int>>> f(n * 2 - 1, vector<vector<int>>(n, vector<int>(n, INT_MIN)));
        f[0][0][0] = grid[0][0];
        for (int k = 1; k < n * 2 - 1; ++k) {
            for (int x1 = max(k - n + 1, 0); x1 <= min(k, n - 1); ++x1) {
                int y1 = k - x1;
                if (grid[x1][y1] == -1) {
                    continue;
                }
                for (int x2 = x1; x2 <= min(k, n - 1); ++x2) {
                    int y2 = k - x2;
                    if (grid[x2][y2] == -1) {
                        continue;
                    }
                    int res = f[k - 1][x1][x2]; // 都往右
                    if (x1) {
                        res = max(res, f[k - 1][x1 - 1][x2]); // 往下，往右
                    }
                    if (x2) {
                        res = max(res, f[k - 1][x1][x2 - 1]); // 往右，往下
                    }
                    if (x1 && x2) {
                        res = max(res, f[k - 1][x1 - 1][x2 - 1]); // 都往下
                    }
                    res += grid[x1][y1];
                    if (x2 != x1) { // 避免重复摘同一个樱桃
                        res += grid[x2][y2];
                    }
                    f[k][x1][x2] = res;
                }
            }
        }
        return max(f.back().back().back(), 0);
    }
}

void cherryPickup_test() {
    vector<vector<int>> grid;
    grid = {{0, 1, -1},
            {1, 0, -1},
            {1, 1, 1}};
    cout << cherryPickup::cherryPickup(grid) << endl;
    grid = {{1,  1,  -1},
            {1,  -1, 1},
            {-1, 1,  1}};
    cout << cherryPickup::cherryPickup(grid) << endl;
}

namespace networkDelayTime {
    int networkDelayTime(vector<vector<int>> &times, int n, int k) {
        const int inf = INT_MAX / 2;
        vector<vector<int>> g(n, vector<int>(n, inf));
        for (auto &t : times) {
            int x = t[0] - 1, y = t[1] - 1;
            g[x][y] = t[2];
        }

        vector<int> dist(n, inf);
        dist[k - 1] = 0;
        vector<int> used(n);
        for (int i = 0; i < n; ++i) {
            int x = -1;
            for (int y = 0; y < n; ++y) {
                if (!used[y] && (x == -1 || dist[y] < dist[x])) {
                    x = y;
                }
            }
            used[x] = true;
            for (int y = 0; y < n; ++y) {
                dist[y] = min(dist[y], dist[x] + g[x][y]);
            }
        }

        int ans = *max_element(dist.begin(), dist.end());
        return ans == inf ? -1 : ans;
    }
}

void networkDelayTime_test() {
    vector<vector<int>> times;
    int n, k;
    times = {
            {2, 1, 1},
            {2, 3, 1},
            {3, 4, 1}
    };
    n = 4, k = 2;
    cout << networkDelayTime::networkDelayTime(times, n, k) << endl;
    times = {{1, 2, 1}};
    n = 2, k = 1;
    cout << networkDelayTime::networkDelayTime(times, n, k) << endl;
}

namespace nextGreatestLetter {
    char nextGreatestLetter(vector<char> &letters, char target) {
        return target < letters.back() ? *upper_bound(letters.begin(), letters.end() - 1, target) : letters[0];
    }
}

void nextGreatestLetter_test() {
    vector<char> letters;
    char target;
    letters = {'c', 'f', 'j'};
    target = 'a';
    cout << nextGreatestLetter::nextGreatestLetter(letters, target) << endl;
    letters = {'c', 'f', 'j'};
    target = 'c';
    cout << nextGreatestLetter::nextGreatestLetter(letters, target) << endl;
    letters = {'x', 'x', 'y', 'y'};
    target = 'z';
    cout << nextGreatestLetter::nextGreatestLetter(letters, target) << endl;
}

namespace minCostClimbingStairs {
    int minCostClimbingStairs(vector<int> &cost) {
        int n = cost.size();
        vector<int> dp0(n), dp1(n);
        bool is_step = false;
        dp0[0] = cost[0];
        dp1[0] = 0;
        for (int i = 1; i < n; ++i) {
            dp0[i] = min(dp0[i - 1], dp1[i - 1]) + cost[i];
            dp1[i] = dp0[i - 1];
        }
        return min(dp0[n - 1], dp1[n - 1]);
    }
}

void minCostClimbingStairs_test() {
    vector<int> cost;
    cost = {10, 15, 20};
    cout << minCostClimbingStairs::minCostClimbingStairs(cost) << endl;
    cost = {1, 100, 1, 1, 1, 100, 1, 1, 100, 1};
    cout << minCostClimbingStairs::minCostClimbingStairs(cost) << endl;
}

namespace dominantIndex {
    int dominantIndex(vector<int> &nums) {
        int index = 0;
        int max = nums[0], secmax = 0;
        for (int i = 1; i < nums.size(); ++i) {
            int num = nums[i];
            if (num > max) {
                secmax = max;
                max = num;
                index = i;
            } else if (num > secmax) {
                secmax = num;
            }
        }
        return secmax * 2 > max ? -1 : index;
    }
}

void dominantIndex_test() {
    vector<int> nums;
    nums = {3, 6, 1, 0};
    cout << dominantIndex::dominantIndex(nums) << endl;
    nums = {1, 2, 3, 4};
    cout << dominantIndex::dominantIndex(nums) << endl;
    nums = {0, 0, 3, 2};
    cout << dominantIndex::dominantIndex(nums) << endl;
}

namespace shortestCompletingWord {
    string shortestCompletingWord(string licensePlate, vector<string> &words) {
        array<int, 26> cnt{};
        for (char ch : licensePlate) {
            if (isalpha(ch)) {
                ++cnt[tolower(ch) - 'a'];
            }
        }
        int idx = -1;
        for (int i = 0; i < words.size(); ++i) {
            array<int, 26> c{};
            for (char ch : words[i]) {
                ++c[ch - 'a'];
            }
            bool ok = true;
            for (int j = 0; j < 26; ++j) {
                if (c[j] < cnt[j]) {
                    ok = false;
                    break;
                }
            }
            if (ok && (idx < 0 || words[i].length() < words[idx].length())) {
                idx = i;
            }
        }
        return words[idx];
    }
}

void shortestCompletingWord_test() {
    string licensePlate;
    vector<string> words;
    licensePlate = "1s3 PSt";
    words = {"step", "steps", "stripe", "stepple"};
    cout << shortestCompletingWord::shortestCompletingWord(licensePlate, words) << endl;
    licensePlate = "1s3 456";
    words = {"looks", "pest", "stew", "show"};
    cout << shortestCompletingWord::shortestCompletingWord(licensePlate, words) << endl;
    licensePlate = "Ah71752";
    words = {"suggest", "letter", "of", "husband", "easy", "education", "drug", "prevent", "writer", "old"};
    cout << shortestCompletingWord::shortestCompletingWord(licensePlate, words) << endl;
}

namespace openLock {
    int openLock(vector<string> &deadends, string target) {
        stack<string> stk;
        if (target == "0000") {
            return 0;
        }
        unordered_set<string> dead(deadends.begin(), deadends.end());
        if (dead.count("0000")) {
            return -1;
        }

        auto num_prev = [](char x) -> char {
            return (x == '0' ? '9' : x - 1);
        };
        auto num_succ = [](char x) -> char {
            return (x == '9' ? '0' : x + 1);
        };

        // 枚举 status 通过一次旋转得到的数字
        auto get = [&](string &status) -> vector<string> {
            vector<string> ret;
            for (int i = 0; i < 4; ++i) {
                char num = status[i];
                status[i] = num_prev(num);
                ret.push_back(status);
                status[i] = num_succ(num);
                ret.push_back(status);
                status[i] = num;
            }
            return ret;
        };
        queue<pair<string, int>> q;
        q.emplace("0000", 0);
        unordered_set<string> seen = {"0000"};

        while (!q.empty()) {
            auto[status, step] = q.front();
            q.pop();
            for (auto &&next_status: get(status)) {
                if (!seen.count(next_status) && !dead.count(next_status)) {
                    if (next_status == target) {
                        return step + 1;
                    }
                    q.emplace(next_status, step + 1);
                    seen.insert(move(next_status));
                }
            }
        }
        return -1;
    }
}

void openLock_test() {
    vector<string> deadends;
    string target;
    deadends = {"0201", "0101", "0102", "1212", "2002"};
    target = "0202";
    cout << openLock::openLock(deadends, target) << endl;
    deadends = {"8888"};
    target = "0009";
    cout << openLock::openLock(deadends, target) << endl;
    deadends = {"8887", "8889", "8878", "8898", "8788", "8988", "7888", "9888"};
    target = "8888";
    cout << openLock::openLock(deadends, target) << endl;
}

namespace crackSafe {
    unordered_set<int> seen;
    string ans;
    int highest;
    int k;

    void dfs(int node) {
        for (int x = 0; x < k; ++x) {
            int nei = node * 10 + x;
            if (!seen.count(nei)) {
                seen.insert(nei);
                dfs(nei % highest);
                ans += (x + '0');
            }
        }
    }

    string crackSafe(int n, int _k) {
        highest = pow(10, n - 1);
        k = _k;
        dfs(0);
        ans += string(n - 1, '0');
        return ans;
    }
}

void crackSafe_test() {
    int n = 1, k = 2;
    cout << crackSafe::crackSafe(n, k) << endl;
    n = 2, k = 2;
    cout << crackSafe::crackSafe(n, k) << endl;
}

namespace reachNumber {
    int reachNumber(int target) {
        target = abs(target);
        int k = 0;
        while (target > 0) {
            k++;
            target -= k;
        }
        return target % 2 == 0 ? k : k + 1 + k % 2;
    }
}

void reachNumber_test() {
    int target;
    target = 2;
    cout << reachNumber::reachNumber(target) << endl;
    target = 3;
    cout << reachNumber::reachNumber(target) << endl;
}

namespace intersectionSizeTwo {
    void help(vector<vector<int>> &intervals, vector<vector<int>> &temp, int pos, int num) {
        for (int i = pos; i >= 0; i--) {
            if (intervals[i][1] < num) {
                break;
            }
            temp[i].push_back(num);
        }
    }

    int intersectionSizeTwo(vector<vector<int>> &intervals) {
        int n = intervals.size();
        int res = 0;
        int m = 2;
        sort(intervals.begin(), intervals.end(), [&](vector<int> &a, vector<int> &b) {
            if (a[0] == b[0]) {
                return a[1] > b[1];
            }
            return a[0] < b[0];
        });
        vector<vector<int>> temp(n);
        for (int i = n - 1; i >= 0; i--) {
            for (int j = intervals[i][0], k = temp[i].size(); k < m; j++, k++) {
                res++;
                help(intervals, temp, i - 1, j);
            }
        }
        return res;
    }
}

void intersectionSizeTwo_test() {
    vector<vector<int>> intervals;
    intervals = {{1, 3},
                 {3, 7},
                 {8, 9}};
    cout << intersectionSizeTwo::intersectionSizeTwo(intervals) << endl;
    intervals = {{1, 3},
                 {1, 4},
                 {2, 5},
                 {3, 5}};
    cout << intersectionSizeTwo::intersectionSizeTwo(intervals) << endl;
    intervals = {{1, 2},
                 {2, 3},
                 {2, 4},
                 {4, 5}};
    cout << intersectionSizeTwo::intersectionSizeTwo(intervals) << endl;
}

namespace makeLargestSpecial {
    string makeLargestSpecial(string s) {
        if (s.size() <= 2) {
            return s;
        }
        int cnt = 0, left = 0;
        vector<string> subs;
        for (int i = 0; i < s.size(); ++i) {
            if (s[i] == '1') {
                ++cnt;
            } else {
                --cnt;
                if (cnt == 0) {
                    subs.push_back("1" + makeLargestSpecial(s.substr(left + 1, i - left - 1)) + "0");
                    left = i + 1;
                }
            }
        }

        sort(subs.begin(), subs.end(), greater<string>{});
        string ans = accumulate(subs.begin(), subs.end(), ""s);
        return ans;
    }
}

void makeLargestSpecial_test() {
    string s, ans;
    s = "11011000";
    ans = makeLargestSpecial::makeLargestSpecial(s);
    cout << ans << endl;
    s = "10";
    ans = makeLargestSpecial::makeLargestSpecial(s);
    cout << ans << endl;
}

namespace countPrimeSetBits {
    int countPrimeSetBits(int left, int right) {
        int ans = 0;
        for (int x = left; x <= right; ++x) {
            if ((1 << __builtin_popcount(x)) & 665772) {
                ++ans;
            }
        }
        return ans;
    }
}


void countPrimeSetBits_test() {
    int left, right;
    left = 6, right = 10;
    cout << countPrimeSetBits::countPrimeSetBits(left, right) << endl;
    left = 10, right = 15;
    cout << countPrimeSetBits::countPrimeSetBits(left, right) << endl;
}

namespace partitionLabels {
    vector<int> partitionLabels(string s) {
        vector<int> last(26);
        int n = s.size();
        for (int i = 0; i < n; ++i) {
            last[s[i] - 'a'] = i;
        }
        vector<int> partition;
        int start = 0, end = 0;
        for (int i = 0; i < n; ++i) {
            end = max(end, last[s[i] - 'a']);
            if (i == end) {
                partition.push_back(end - start + 1);
                start = end + 1;
            }
        }
        return partition;
    }
}

void partitionLabels_test() {
    string s;
    vector<int> ans;
    s = "ababcbacadefegdehijhklij";
    ans = partitionLabels::partitionLabels(s);
    print_vector(ans);
    s = "eccbbbbdec";
    ans = partitionLabels::partitionLabels(s);
    print_vector(ans);
    s = "caedbdedda";
    ans = partitionLabels::partitionLabels(s);
    print_vector(ans);
    s = "eaaaabaaec";
    ans = partitionLabels::partitionLabels(s);
    print_vector(ans);
}

namespace orderOfLargestPlusSign {
    int orderOfLargestPlusSign(int n, vector<vector<int>> &mines) {
        vector<vector<int>> dp(n, vector<int>(n, n));
        unordered_set<int> banned;
        for (auto &&vec : mines) {
            banned.emplace(vec[0] * n + vec[1]);
        }
        int ans = 0;
        for (int i = 0; i < n; i++) {
            int count = 0;
            /* left */
            for (int j = 0; j < n; j++) {
                if (banned.count(i * n + j)) {
                    count = 0;
                } else {
                    count++;
                }
                dp[i][j] = min(dp[i][j], count);
            }
            count = 0;
            /* right */
            for (int j = n - 1; j >= 0; j--) {
                if (banned.count(i * n + j)) {
                    count = 0;
                } else {
                    count++;
                }
                dp[i][j] = min(dp[i][j], count);
            }
        }
        for (int i = 0; i < n; i++) {
            int count = 0;
            /* up */
            for (int j = 0; j < n; j++) {
                if (banned.count(j * n + i)) {
                    count = 0;
                } else {
                    count++;
                }
                dp[j][i] = min(dp[j][i], count);
            }
            count = 0;
            /* down */
            for (int j = n - 1; j >= 0; j--) {
                if (banned.count(j * n + i)) {
                    count = 0;
                } else {
                    count++;
                }
                dp[j][i] = min(dp[j][i], count);
                ans = max(ans, dp[j][i]);
            }
        }
        return ans;
    }
}

void orderOfLargestPlusSign_test() {
    int n;
    vector<vector<int>> mines;
    n = 5, mines = {{4, 2}};
    cout << orderOfLargestPlusSign::orderOfLargestPlusSign(n, mines) << endl;
    n = 1, mines = {{0, 0}};
    cout << orderOfLargestPlusSign::orderOfLargestPlusSign(n, mines) << endl;
    n = 2, mines = {{0, 0},
                    {0, 1},
                    {1, 0}};
    cout << orderOfLargestPlusSign::orderOfLargestPlusSign(n, mines) << endl;
    n = 3, mines = {{0, 1},
                    {0, 2},
                    {1, 0},
                    {1, 1},
                    {1, 2},
                    {2, 0},
                    {2, 1},
                    {2, 2}};
    cout << orderOfLargestPlusSign::orderOfLargestPlusSign(n, mines) << endl;
}

namespace minSwapsCouples {
    int minSwapsCouples(vector<int> &row) {
        int n = row.size();
        int tot = n / 2;

        vector<vector<int>> graph(tot);
        for (int i = 0; i < n; i += 2) {
            int l = row[i] / 2;
            int r = row[i + 1] / 2;
            if (l != r) {
                graph[l].push_back(r);
                graph[r].push_back(l);
            }
        }
        vector<int> visited(tot, 0);
        int ret = 0;
        for (int i = 0; i < tot; i++) {
            if (visited[i] == 0) {
                queue<int> q;
                visited[i] = 1;
                q.push(i);
                int cnt = 0;

                while (!q.empty()) {
                    int x = q.front();
                    q.pop();
                    cnt += 1;

                    for (int y: graph[x]) {
                        if (visited[y] == 0) {
                            visited[y] = 1;
                            q.push(y);
                        }
                    }
                }
                ret += cnt - 1;
            }
        }
        return ret;
    }
}

void minSwapsCouples_test() {
    vector<int> row;
    row = {0, 2, 1, 3};
    cout << minSwapsCouples::minSwapsCouples(row) << endl;
    row = {3, 2, 0, 1};
    cout << minSwapsCouples::minSwapsCouples(row) << endl;
}

namespace isToeplitzMatrix {
    bool isToeplitzMatrix(vector<vector<int>> &matrix) {
        int m = matrix.size(), n = matrix[0].size();
        for (int i = 1; i < m; i++) {
            for (int j = 1; j < n; j++) {
                if (matrix[i][j] != matrix[i - 1][j - 1]) {
                    return false;
                }
            }
        }
        return true;
    }
}

void isToeplitzMatrix_test() {
    vector<vector<int>> matrix;
    matrix = {{1, 2, 3, 4},
              {5, 1, 2, 3},
              {9, 5, 1, 2}};
    cout << isToeplitzMatrix::isToeplitzMatrix(matrix) << endl;
    matrix = {{1, 2},
              {2, 2}};
    cout << isToeplitzMatrix::isToeplitzMatrix(matrix) << endl;
}

namespace reorganizeString {
    string reorganizeString(string s) {
        if (s.length() < 2) {
            return s;
        }
        vector<int> counts(26, 0);
        int maxCount = 0;
        int length = s.length();
        for (int i = 0; i < length; i++) {
            char c = s[i];
            counts[c - 'a']++;
            maxCount = max(maxCount, counts[c - 'a']);
        }
        if (maxCount > (length + 1) / 2) {
            return "";
        }
        string reorganizeArray(length, ' ');
        int evenIndex = 0, oddIndex = 1;
        int halfLength = length / 2;
        for (int i = 0; i < 26; i++) {
            char c = 'a' + i;
            while (counts[i] > 0 && counts[i] <= halfLength && oddIndex < length) {
                reorganizeArray[oddIndex] = c;
                counts[i]--;
                oddIndex += 2;
            }
            while (counts[i] > 0) {
                reorganizeArray[evenIndex] = c;
                counts[i]--;
                evenIndex += 2;
            }
        }
        return reorganizeArray;
    }
}

void reorganizeString_test() {
    string s;
    s = "baa";
    cout << reorganizeString::reorganizeString(s) << endl;
    cout << "--------" << endl;
    s = "aab";
    cout << reorganizeString::reorganizeString(s) << endl;
    cout << "--------" << endl;
    s = "aaab";
    cout << reorganizeString::reorganizeString(s) << endl;
    cout << "--------" << endl;
}

namespace maxChunksToSorted {
    int maxChunksToSorted(vector<int> &arr) {
        unordered_map<int, int> cnt;
        int res = 0;
        vector<int> sortedArr = arr;
        sort(sortedArr.begin(), sortedArr.end());
        for (int i = 0; i < sortedArr.size(); i++) {
            int x = arr[i], y = sortedArr[i];
            cnt[x]++;
            if (cnt[x] == 0) {
                cnt.erase(x);
            }
            cnt[y]--;
            if (cnt[y] == 0) {
                cnt.erase(y);
            }
            if (cnt.size() == 0) {
                res++;
            }
        }
        return res;
    }
};

void maxChunksToSorted_test() {
    vector<int> arr;
    arr = {5, 4, 3, 2, 1};
    cout << maxChunksToSorted::maxChunksToSorted(arr) << endl;
    arr = {2, 1, 3, 4, 4};
    cout << maxChunksToSorted::maxChunksToSorted(arr) << endl;
}

namespace numJewelsInStones {
    int numJewelsInStones(string jewels, string stones) {
        unordered_set<char> jewels_set;
        int ans = 0;
        for (auto c : jewels) {
            jewels_set.insert(c);
        }
        for (auto c : stones) {
            if (jewels_set.count(c)) {
                ans++;
            }
        }
        return ans;
    }
}

void numJewelsInStones_test() {
    string jewels = "aA", stones = "aAAbbbb";
    cout << numJewelsInStones::numJewelsInStones(jewels, stones) << endl;
    jewels = "Z", stones = "zz";
    cout << numJewelsInStones::numJewelsInStones(jewels, stones) << endl;
}

namespace slidingPuzzle {
    vector<vector<int>> neighbors = {{1, 3},
                                     {0, 2, 4},
                                     {1, 5},
                                     {0, 4},
                                     {1, 3, 5},
                                     {2, 4}};

    int slidingPuzzle(vector<vector<int>> &board) {
        // 枚举status，通过一次交换操作得到的状态
        auto get = [&](string &status) -> vector<string> {
            vector<string> ret;
            int x = status.find('0');
            for (int y: neighbors[x]) {
                swap(status[x], status[y]);
                ret.push_back(status);
                swap(status[x], status[y]);
            }
            return ret;
        };
        string initial;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 3; ++j) {
                initial += char(board[i][j] + '0');
            }
        }
        if (initial == "123450") {
            return 0;
        }

        queue<pair<string, int>> q;
        q.emplace(initial, 0);
        unordered_set<string> seen = {initial};

        while (!q.empty()) {
            auto[status, step] = q.front();
            q.pop();
            for (auto &&next_status: get(status)) {
                if (!seen.count(next_status)) {
                    if (next_status == "123450") {
                        return step + 1;
                    }
                    q.emplace(next_status, step + 1);
                    seen.insert(move(next_status));
                }
            }
        }

        return -1;
    }
}

void slidingPuzzle_test() {
    vector<vector<int>> board;
    board = {{1, 2, 3},
             {4, 0, 5}};
    cout << slidingPuzzle::slidingPuzzle(board) << endl;
    board = {{1, 2, 3},
             {5, 4, 0}};
    cout << slidingPuzzle::slidingPuzzle(board) << endl;
    board = {{4, 1, 2},
             {5, 0, 3}};
    cout << slidingPuzzle::slidingPuzzle(board) << endl;
}

namespace isIdealPermutation {
    bool isIdealPermutation(vector<int> &nums) {
        int n = nums.size(), minSuff = nums[n - 1];
        for (int i = n - 3; i >= 0; i--) {
            if (nums[i] > minSuff) {
                return false;
            }
            minSuff = min(minSuff, nums[i + 1]);
        }
        return true;

    }
}

void isIdealPermutation_test() {
    vector<int> nums;
    nums = {1, 0, 2};
    cout << isIdealPermutation::isIdealPermutation(nums) << endl;
    nums = {1, 2, 0};
    cout << isIdealPermutation::isIdealPermutation(nums) << endl;
}

namespace canTransform {
    bool canTransform(string start, string end) {
        int n = start.length();
        int i = 0, j = 0;
        while (i < n && j < n) {
            while (i < n && start[i] == 'X') {
                i++;
            }
            while (j < n && end[j] == 'X') {
                j++;
            }
            if (i < n && j < n) {
                if (start[i] != end[j]) {
                    return false;
                }
                char c = start[i];
                if ((c == 'L' && i < j) || (c == 'R' && i > j)) {
                    return false;
                }
                i++;
                j++;
            }
        }
        while (i < n) {
            if (start[i] != 'X') {
                return false;
            }
            i++;
        }
        while (j < n) {
            if (end[j] != 'X') {
                return false;
            }
            j++;
        }
        return true;
    }
}

void canTransform_test() {
    string start, end;
    start = "RXXLRXRXL";
    end = "XRLXXRRLX";
    cout << canTransform::canTransform(start, end) << endl;
    start = "X";
    end = "L";
    cout << canTransform::canTransform(start, end) << endl;
}

namespace kthGrammar {
    int kthGrammar(int n, int k) {
        if (n == 1)
            return 0;
        return (k & 1) ^ 1 ^ kthGrammar(n - 1, (k + 1) / 2);
    }
}

void kthGrammar_test() {
    int n, k;
    n = 1, k = 1;
    cout << kthGrammar::kthGrammar(n, k) << endl;
    n = 2, k = 1;
    cout << kthGrammar::kthGrammar(n, k) << endl;
    n = 2, k = 2;
    cout << kthGrammar::kthGrammar(n, k) << endl;
}

namespace reachingPoints {
    bool reachingPoints(int sx, int sy, int tx, int ty) {
        while (tx > sx && ty > sy && tx != ty) {
            if (tx > ty) {
                tx %= ty;
            } else {
                ty %= tx;
            }
        }

        if (tx == sx && ty == sy) {
            return true;
        } else if (tx == sx) {
            return ty > sy && (ty - sy) % tx == 0;
        } else if (ty == sy) {
            return tx > sx && (tx - sx) % ty == 0;
        } else
            return false;
        return false;
    }
}

void reachingPoints_test() {
    int sx, sy, tx, ty;
    sx = 1, sy = 1, tx = 3, ty = 5;
    cout << reachingPoints::reachingPoints(sx, sy, tx, ty) << endl;
    sx = 1, sy = 1, tx = 2, ty = 2;
    cout << reachingPoints::reachingPoints(sx, sy, tx, ty) << endl;
    sx = 1, sy = 1, tx = 1, ty = 1;
    cout << reachingPoints::reachingPoints(sx, sy, tx, ty) << endl;
}

namespace numRabbits {
    int numRabbits(vector<int> &answers) {
        unordered_map<int, int> map;
        for (auto ans:answers) {
            ++map[ans];
        }
        int ans = 0;
        for (auto &[y, x]: map) {
            ans += (x + y) / (y + 1) * (y + 1);
        }
        return ans;
    }
}

void numRabbits_test() {
    vector<int> answers;
    answers = {1, 1, 2};
    cout << numRabbits::numRabbits(answers) << endl;
    answers = {10, 10, 10};
    cout << numRabbits::numRabbits(answers) << endl;
}

namespace minDiffInBST {
    void dfs(TreeNode::TreeNode *root, int &pre, int &ans) {
        if (root == nullptr)
            return;
        dfs(root->left, pre, ans);
        if (pre == -1) {
            pre = root->val;
        } else {
            int tmp = fabs(root->val - pre);
            ans = min(ans, tmp);
            pre = root->val;
        }
        dfs(root->right, pre, ans);
    }

    int minDiffInBST(TreeNode::TreeNode *root) {
        int ans = INT_MAX, pre = -1;
        dfs(root, pre, ans);
        return ans;
    }
}

void minDiffInBST_test() {
    TreeNode::TreeNode *root;
    vector<int> tree;
    tree = {4, 2, 6, 1, 3};
    root = create_treenode(tree, true);
    cout << minDiffInBST::minDiffInBST(root) << endl;
    tree = {1, 0, 48, -1, -1, 12, 4};
    root = create_treenode(tree, true);
    cout << minDiffInBST::minDiffInBST(root) << endl;
}

namespace letterCasePermutation {
    vector<string> letterCasePermutation(string s) {
        vector<string> ans;
        queue<string> qu;
        qu.emplace("");
        while (!qu.empty()) {
            string &curr = qu.front();
            if (curr.size() == s.size()) {
                ans.emplace_back(curr);
                qu.pop();
            } else {
                int pos = curr.size();
                if (isalpha(s[pos])) {
                    string next = curr;
                    next.push_back(s[pos] ^ 32);
                    qu.emplace(next);
                }
                curr.push_back(s[pos]);
            }
        }
        return ans;
    }
}

void letterCasePermutation_test() {
    string s;
    vector<string> ans;
    s = "a1b2";
    ans = letterCasePermutation::letterCasePermutation(s);
    print_vector(ans);
    s = "3z4";
    ans = letterCasePermutation::letterCasePermutation(s);
    print_vector(ans);
}

namespace isBipartite {
    static constexpr int UNCOLORED = 0;
    static constexpr int RED = 1;
    static constexpr int GREEN = 2;
    vector<int> color;
    bool valid;

    void dfs(int node, int c, const vector<vector<int>> &graph) {
        color[node] = c;
        int cNei = (c == RED ? GREEN : RED);
        for (int neighbor: graph[node]) {
            if (color[neighbor] == UNCOLORED) {
                dfs(neighbor, cNei, graph);
                if (!valid) {
                    return;
                }
            } else if (color[neighbor] != cNei) {
                valid = false;
                return;
            }
        }
    }

    bool isBipartite(vector<vector<int>> &graph) {
        int n = graph.size();
        valid = true;
        color.assign(n, UNCOLORED);
        for (int i = 0; i < n && valid; ++i) {
            if (color[i] == UNCOLORED) {
                dfs(i, RED, graph);
            }
        }
        return valid;
    }
}

void isBipartite_test() {
    vector<vector<int>> graph;
    graph = {{1, 2, 3},
             {0, 2},
             {0, 1, 3},
             {0, 2}};
    cout << isBipartite::isBipartite(graph) << endl;
    graph = {{1, 3},
             {0, 2},
             {1, 3},
             {0, 2}};
    cout << isBipartite::isBipartite(graph) << endl;
}

namespace kthSmallestPrimeFraction {
    vector<int> kthSmallestPrimeFraction(vector<int> &arr, int k) {
        int n = arr.size();
        vector<pair<int, int>> frac;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                frac.emplace_back(arr[i], arr[j]);
            }
        }
        sort(frac.begin(), frac.end(), [&](const auto &x, const auto &y) {
            return x.first * y.second < x.second * y.first;
        });
        return {frac[k - 1].first, frac[k - 1].second};
    }
}

void kthSmallestPrimeFraction_test() {
    vector<int> arr, ans;
    int k;
    arr = {1, 2, 3, 5};
    k = 3;
    ans = kthSmallestPrimeFraction::kthSmallestPrimeFraction(arr, k);
    print_vector(ans);
    arr = {1, 7};
    k = 1;
    ans = kthSmallestPrimeFraction::kthSmallestPrimeFraction(arr, k);
    print_vector(ans);
}

namespace findCheapestPrice {
    static constexpr int INF = 10000 * 101 + 1;

    int findCheapestPrice(int n, vector<vector<int>> &flights, int src, int dst, int k) {
        vector<vector<int>> f(k + 2, vector<int>(n, INF));
        f[0][src] = 0;
        for (int t = 1; t <= k + 1; ++t) {
            for (auto &&flight : flights) {
                int j = flight[0], i = flight[1], cost = flight[2];
                f[t][i] = min(f[t][i], f[t - 1][j] + cost);
            }
        }
        int ans = INF;
        for (int t = 1; t <= k + 1; ++t) {
            ans = min(ans, f[t][dst]);
        }
        return (ans == INF ? -1 : ans);
    }
}

void findCheapestPrice_test() {
    int n, k, src, dst;
    vector<vector<int>> flights;
    n = 4, src = 0, dst = 3, k = 1;
    flights = {{0, 1, 100},
               {1, 2, 100},
               {2, 0, 100},
               {1, 3, 600},
               {2, 3, 200}};
    cout << findCheapestPrice::findCheapestPrice(n, flights, src, dst, k) << endl;
    n = 3, src = 0, dst = 2, k = 1;
    flights = {{0, 1, 100},
               {1, 2, 100},
               {0, 2, 500}};
    cout << findCheapestPrice::findCheapestPrice(n, flights, src, dst, k) << endl;
}

namespace rotatedDigits {
    static constexpr int check[10] = {0, 0, 1, -1, -1, 1, 1, -1, 0, 1};

    int rotatedDigits(int n) {
        int ans = 0;
        for (int i = 1; i <= n; ++i) {
            string num = to_string(i);
            bool valid = true, diff = false;
            for (char ch: num) {
                if (check[ch - '0'] == -1) {
                    valid = false;
                } else if (check[ch - '0'] == 1) {
                    diff = true;
                }
            }
            if (valid && diff) {
                ++ans;
            }
        }
        return ans;
    }
}

void rotatedDigits_test() {
    int n = 10;
    cout << rotatedDigits::rotatedDigits(n) << endl;
}

namespace escapeGhosts {
    int manhattanDistance(vector<int> &point1, vector<int> &point2) {
        return abs(point1[0] - point2[0]) + abs(point1[1] - point2[1]);
    }

    bool escapeGhosts(vector<vector<int>> &ghosts, vector<int> &target) {
        vector<int> src(2);
        auto distance = manhattanDistance(src, target);
        for (auto ghost : ghosts) {
            if (manhattanDistance(ghost, target) <= distance) {
                return false;
            }
        }
        return true;
    }
}

void escapeGhosts_test() {
    vector<vector<int>> ghosts;
    vector<int> target;
    ghosts = {{1, 0},
              {0, 3}};
    target = {0, 1};
    cout << escapeGhosts::escapeGhosts(ghosts, target) << endl;
    ghosts = {{1, 0}};
    target = {2, 0};
    cout << escapeGhosts::escapeGhosts(ghosts, target) << endl;
    ghosts = {{2, 0}};
    target = {1, 0};
    cout << escapeGhosts::escapeGhosts(ghosts, target) << endl;
}

namespace numTilings {
    const long long mod = 1e9 + 7;

    int numTilings(int n) {
        vector<vector<long long>> dp(n + 1, vector<long long>(4));
        dp[0][3] = 1;
        for (int i = 1; i <= n; i++) {
            dp[i][0] = dp[i - 1][3];
            dp[i][1] = (dp[i - 1][0] + dp[i - 1][2]) % mod;
            dp[i][2] = (dp[i - 1][0] + dp[i - 1][1]) % mod;
            dp[i][3] = (dp[i - 1][0] + dp[i - 1][1] + dp[i - 1][2] + dp[i - 1][3]) % mod;
        }
        return dp[n][3];
    }
}

void numTilings_test() {
    int n = 3;
    cout << numTilings::numTilings(n) << endl;
    n = 1;
    cout << numTilings::numTilings(n) << endl;
}

namespace customSortString {
    string customSortString(string order, string s) {
        vector<int> val(26);
        for (int i = 0; i < order.size(); ++i) {
            val[order[i] - 'a'] = i + 1;
        }
        sort(s.begin(), s.end(), [&](char c0, char c1) {
            return val[c0 - 'a'] < val[c1 - 'a'];
        });
        return s;
    }
}

void customSortString_test() {
    string order = "cba";
    string s = "abcd";
    cout << customSortString::customSortString(order, s) << endl;
    order = "cbafg";
    s = "abcd";
    cout << customSortString::customSortString(order, s) << endl;
}

namespace numMatchingSubseq {
    int numMatchingSubseq(string s, vector<string> &words) {
        vector<queue<pair<int, int>>> queues(26);
        for (int i = 0; i < words.size(); ++i) {
            queues[words[i][0] - 'a'].emplace(i, 0);
        }
        int res = 0;
        for (char c : s) {
            auto &q = queues[c - 'a'];
            int size = q.size();
            while (size--) {
                auto[i, j] = q.front();
                q.pop();
                ++j;
                if (j == words[i].size()) {
                    ++res;
                } else {
                    queues[words[i][j] - 'a'].emplace(i, j);
                }
            }
        }
        return res;
    }
}

namespace preimageSizeFZF {
    int zeta(long x) {
        int res = 0;
        while (x) {
            res += x / 5;
            x /= 5;
        }
        return res;
    }

    int help(int k) {
        long long r = 5LL * k;
        long long l = 0;
        while (l <= r) {
            long long mid = (l + r) / 2;
            if (zeta(mid) < k) {
                l = mid + 1;
            } else {
                r = mid - 1;
            }
        }
        return r + 1;
    }

    int preimageSizeFZF(int k) {
        return help(k + 1) - help(k);
    }
}

void preimageSizeFZF_test() {
    int k;
    k = 0;
    cout << preimageSizeFZF::preimageSizeFZF(k) << endl;
    k = 5;
    cout << preimageSizeFZF::preimageSizeFZF(k) << endl;
    k = 3;
    cout << preimageSizeFZF::preimageSizeFZF(k) << endl;
}

namespace validTicTacToe {
    bool win(vector<string> board, char p) {
        for (int i = 0; i < 3; ++i) {
            if ((p == board[0][i] && p == board[1][i] && p == board[2][i]) ||
                (p == board[i][0] && p == board[i][1] && p == board[i][2])) {
                return true;
            }
        }
        return ((p == board[0][0] && p == board[1][1] && p == board[2][2]) ||
                (p == board[0][2] && p == board[1][1] && p == board[2][0]));
    }

    bool validTicTacToe(vector<string> &board) {
        int xCount = 0, oCount = 0;
        for (string &row : board) {
            for (char c : row) {
                xCount = (c == 'X') ? (xCount + 1) : xCount;
                oCount = (c == 'O') ? (oCount + 1) : oCount;
            }
        }
        return !((oCount != xCount && oCount != xCount - 1) ||
                 (oCount != xCount - 1 && win(board, 'X')) ||
                 (oCount != xCount && win(board, 'O')));
    }
}

void validTicTacToe_test() {
    vector<string> board;
    board = {"O  ", "   ", "   "};
    cout << validTicTacToe::validTicTacToe(board) << endl;
    board = {"XOX", " X ", "   "};
    cout << validTicTacToe::validTicTacToe(board) << endl;
    board = {"XOX", "O O", "XOX"};
    cout << validTicTacToe::validTicTacToe(board) << endl;
}

namespace numSubarrayBoundedMax {
    int numSubarrayBoundedMax(vector<int> &nums, int left, int right) {
        int res = 0, last2 = -1, last1 = -1;
        for (int i = 0; i < nums.size(); i++) {
            if (nums[i] >= left && nums[i] <= right) {
                last1 = i;
            } else if (nums[i] > right) {
                last2 = i;
                last1 = -1;
            }
            if (last1 != -1) {
                res += last1 - last2;
            }
        }
        return res;
    }
}

void numMatchingSubseq_test() {
    string s;
    vector<string> words;
    s = "abcde";
    words = {"a", "bb", "acd", "ace"};
    cout << numMatchingSubseq::numMatchingSubseq(s, words) << endl;
    s = "dsahjpjauf";
    words = {"ahjpjau", "ja", "ahbwzgqnuk", "tnmlanowax"};
    cout << numMatchingSubseq::numMatchingSubseq(s, words) << endl;
}

void numSubarrayBoundedMax_test() {
    vector<int> nums;
    int left, right;
    nums = {2, 1, 4, 3};
    left = 2, right = 3;
    cout << numSubarrayBoundedMax::numSubarrayBoundedMax(nums, left, right) << endl;
    nums = {2, 9, 2, 5, 6};
    left = 2, right = 8;
    cout << numSubarrayBoundedMax::numSubarrayBoundedMax(nums, left, right) << endl;

}

namespace rotateString {
    bool rotateString(string s, string goal) {
        return s.size() == goal.size() && (s + s).find(goal) != string::npos;
    }
}

void rotateString_test() {
    string s = "abcde";
    string goal = "cdeab";
    cout << rotateString::rotateString(s, goal) << endl;
    s = "abcde", goal = "abced";
    cout << rotateString::rotateString(s, goal) << endl;
}

namespace allPathsSourceTarget {
    vector<vector<int>> ans;
    vector<int> stk;

    void dfs(vector<vector<int>> &graph, int x, int n) {
        if (x == n) {
            ans.push_back(stk);
            return;
        }
        for (auto &y : graph[x]) {
            stk.push_back(y);
            dfs(graph, y, n);
            stk.pop_back();
        }
    }

    vector<vector<int>> allPathsSourceTarget(vector<vector<int>> &graph) {
        ans.clear();
        stk.clear();
        stk.push_back(0);
        dfs(graph, 0, graph.size() - 1);
        return ans;
    }
}

void allPathsSourceTarget_test() {
    vector<vector<int>> paths;
    vector<vector<int>> graph = {{1, 2},
                                 {3},
                                 {3},
                                 {}};
    paths = allPathsSourceTarget::allPathsSourceTarget(graph);
    for (auto path : paths) {
        print_vector(path);
    }
    cout << "++++++++++++++++++" << endl;
    graph = {{4, 3, 1},
             {3, 2, 4},
             {3},
             {4},
             {}};
    paths = allPathsSourceTarget::allPathsSourceTarget(graph);
    for (auto path : paths) {
        print_vector(path);
    }
    cout << "++++++++++++++++++" << endl;
}

namespace bestRotation {
    int bestRotation(vector<int> &nums) {
        int n = nums.size();
        vector<int> diffs(n);
        for (int i = 0; i < n; ++i) {
            int low = (i + 1) % n;
            int high = (i - nums[i] + n + 1) % n;
            diffs[low]++;
            diffs[high]--;
            if (low >= high) {
                diffs[0]++;
            }
        }
        int bestIndex = 0;
        int maxScore = 0;
        int score = 0;
        for (int i = 0; i < n; ++i) {
            score += diffs[i];
            if (score > maxScore) {
                bestIndex = i;
                maxScore = score;
            }
        }
        return bestIndex;
    }
};

void bestRotation_test() {
    vector<int> nums;
    nums = {2, 3, 1, 4, 0};
    cout << bestRotation::bestRotation(nums) << endl;
    nums = {1, 3, 0, 2, 4};
    cout << bestRotation::bestRotation(nums) << endl;

}

namespace champagneTower {
    double champagneTower(int poured, int query_row, int query_glass) {
        vector<double> row = {(double) poured};
        for (int i = 1; i <= query_row; ++i) {
            vector<double> nextRow(i + 1, 0.0);
            for (int j = 0; j < row.size(); ++j) {
                double volume = row[j];
                if (volume > 1) {
                    nextRow[j] += (volume - 1) / 2;
                    nextRow[j + 1] += (volume - 1) / 2;
                }
            }
            row = nextRow;
        }
        return min(1.0, row[query_glass]);
    }
}

void champagneTower_test() {
    int poured, query_row, query_glass;
    poured = 1, query_row = 1, query_glass = 1;
    cout << champagneTower::champagneTower(poured, query_row, query_glass) << endl;
    poured = 2, query_row = 1, query_glass = 1;
    cout << champagneTower::champagneTower(poured, query_row, query_glass) << endl;
    poured = 100000009, query_row = 33, query_glass = 17;
    cout << champagneTower::champagneTower(poured, query_row, query_glass) << endl;
}

namespace minSwap {
    int minSwap(vector<int> &nums1, vector<int> &nums2) {
        int n = nums1.size();
        int a = 0, b = 1;
        for (int i = 1; i < n; ++i) {
            int at = a, bt = b;
            a = b = n;
            if (nums1[i] > nums1[i - 1] && nums2[i] > nums2[i - 1]) {
                a = min(a, at);
                b = min(b, bt + 1);
            }
            if (nums1[i] > nums2[i - 1] && nums2[i] > nums1[i - 1]) {
                a = min(a, bt);
                b = min(b, at + 1);
            }
        }
        return min(a, b);
    }
}

void minSwap_test() {
    vector<int> nums1, nums2;
    nums1 = {1, 3, 5, 4};
    nums2 = {1, 2, 3, 7};
    cout << minSwap::minSwap(nums1, nums2) << endl;
    nums1 = {0, 3, 5, 8, 9};
    nums2 = {2, 1, 4, 6, 9};
    cout << minSwap::minSwap(nums1, nums2) << endl;
}

namespace uniqueMorseRepresentations {
    int uniqueMorseRepresentations(vector<string> &words) {
        vector<string> code_table;
        code_table = {".-", "-...", "-.-.", "-..", ".", "..-.", "--.", "....", "..", ".---", "-.-", ".-..", "--", "-.",
                      "---", ".--.", "--.-", ".-.", "...", "-", "..-", "...-", ".--", "-..-", "-.--", "--.."};
        unordered_map<char, string> map;
        for (int i = 0; i < code_table.size(); ++i) {
            map[i + 'a'] = code_table[i];
        }
        unordered_set<string> code_set;
        for (auto word : words) {
            string s;
            for (auto c : word) {
                s = s + map[c];
            }
            code_set.insert(s);
        }
        return code_set.size();
    }
}

void uniqueMorseRepresentations_test() {
    vector<string> words;
    words = {"gin", "zen", "gig", "msg"};
    cout << uniqueMorseRepresentations::uniqueMorseRepresentations(words) << endl;
    words = {"a"};
    cout << uniqueMorseRepresentations::uniqueMorseRepresentations(words) << endl;
}

namespace splitArraySameAverage {
    bool splitArraySameAverage(vector<int> &nums) {
        int n = nums.size(), m = n / 2;
        int sum = accumulate(nums.begin(), nums.end(), 0);
        bool isPossible = false;
        for (int i = 1; i <= m; ++i) {
            if (sum * i % n == 0) {
                isPossible = true;
                break;
            }
        }
        if (!isPossible) {
            return false;
        }
        vector<unordered_set<int>> dp(m + 1);
        dp[0].insert(0);
        for (int num : nums) {
            for (int i = m; i >= 1; --i) {
                for (int x : dp[i - 1]) {
                    int curr = x + num;
                    if (curr * n == sum * i) {
                        return true;
                    }
                    dp[i].emplace(curr);
                }
            }
        }
        return false;
    }
}

void splitArraySameAverage_test() {
    vector<int> nums;
    nums = {1, 2, 3, 4, 5, 6, 7, 8};
    cout << splitArraySameAverage::splitArraySameAverage(nums) << endl;
    nums = {3, 1};
    cout << splitArraySameAverage::splitArraySameAverage(nums) << endl;
}

namespace numberOfLines {
    vector<int> numberOfLines(vector<int> &widths, string s) {
        int line = 100;
        int cur = 0;
        int ansline = 0;
        for (auto c : s) {
            if (cur + widths[c - 'a'] <= 100) {
                cur += widths[c - 'a'];
            } else {
                ansline++;
                cur = widths[c - 'a'];
            }
        }
        ansline = cur > 0 ? ansline + 1 : ansline;
        return {ansline, cur};
    }
}

void numberOfLines_test() {
    vector<int> widths, ans;
    string s;
    s = "abcdefghijklmnopqrstuvwxyz";
    widths = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    ans = numberOfLines::numberOfLines(widths, s);
    print_vector(ans);
    cout << "-------------" << endl;
    widths = {4, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    s = "bbbcccdddaaa";
    ans = numberOfLines::numberOfLines(widths, s);
    print_vector(ans);
    cout << "-------------" << endl;
}

namespace maxIncreaseKeepingSkyline {
    int maxIncreaseKeepingSkyline(vector<vector<int>> &grid) {
        int n = grid.size();
        vector<int> rowMax(n);
        vector<int> colMax(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                rowMax[i] = max(rowMax[i], grid[i][j]);
                colMax[j] = max(colMax[j], grid[i][j]);
            }
        }
        int ans = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                ans += min(rowMax[i], colMax[j]) - grid[i][j];
            }
        }
        return ans;
    }
}

void maxIncreaseKeepingSkyline_test() {
    vector<vector<int>> grid;
    grid = {{3, 0, 8, 4},
            {2, 4, 5, 7},
            {9, 2, 6, 3},
            {0, 3, 1, 0}};
    cout << maxIncreaseKeepingSkyline::maxIncreaseKeepingSkyline(grid) << endl;
    grid = {{0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};
    cout << maxIncreaseKeepingSkyline::maxIncreaseKeepingSkyline(grid) << endl;
}

namespace xorGame {
    bool xorGame(vector<int> &nums) {
        if (nums.size() % 2 == 0) {
            return true;
        }
        int xorsum = 0;
        for (int num : nums) {
            xorsum ^= num;
        }
        return xorsum == 0;
    }
}

void xorGame_test() {
    vector<int> nums;
    nums = {1, 1, 2};
    cout << xorGame::xorGame(nums) << endl;
    nums = {1, 2};
    cout << xorGame::xorGame(nums) << endl;
    nums = {1, 2, 3};
    cout << xorGame::xorGame(nums) << endl;
}

namespace largestTriangleArea {
    double triangleArea(int x1, int y1, int x2, int y2, int x3, int y3) {
        return 0.5 * abs(x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2);
    }

    double largestTriangleArea(vector<vector<int>> &points) {
        int n = points.size();
        double ret = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++) {
                    ret = max(ret, triangleArea(points[i][0], points[i][1], points[j][0], points[j][1], points[k][0],
                                                points[k][1]));
                }
            }
        }
        return ret;
    }
}

void largestTriangleArea_test() {
    vector<vector<int>> points;
    points = {{0, 1},
              {1, 0},
              {0, 2},
              {2, 0}};
    cout << largestTriangleArea::largestTriangleArea(points) << endl;
    points = {{1, 0},
              {0, 0},
              {0, 1}};
    cout << largestTriangleArea::largestTriangleArea(points) << endl;
    points = {{4, 6},
              {6, 5},
              {3, 1}};
    cout << largestTriangleArea::largestTriangleArea(points) << endl;
};

namespace largestSumOfAverages {
    double largestSumOfAverages(vector<int> &nums, int k) {
        int n = nums.size();
        vector<double> prefix(n + 1);
        for (int i = 0; i < n; ++i) {
            prefix[i + 1] = prefix[i] + nums[i];
        }
        vector<vector<double>> dp(n + 1, vector<double>(k + 1));
        for (int i = 1; i <= n; i++) {
            dp[i][1] = prefix[i] / i;
        }
        for (int j = 2; j <= k; ++j) {
            for (int i = j; i <= n; ++i) {
                for (int x = j - 1; x < i; ++x) {
                    dp[i][j] = max(dp[i][j], dp[x][j - 1] + (prefix[i] - prefix[x]) / (i - x));
                }
            }
        }
        return dp[n][k];
    }
}

void largestSumOfAverages_test() {
    vector<int> nums;
    int k;
    k = 3;
    nums = {9, 1, 2, 3, 9};
    cout << largestSumOfAverages::largestSumOfAverages(nums, k) << endl;
    nums = {1, 2, 3, 4, 5, 6, 7};
    k = 4;
    cout << largestSumOfAverages::largestSumOfAverages(nums, k) << endl;
}

namespace pruneTree {
    TreeNode::TreeNode *pruneTree(TreeNode::TreeNode *root) {
        if (!root) {
            return nullptr;
        }
        root->left = pruneTree(root->left);
        root->right = pruneTree(root->right);
        if (!root->left && !root->right && !root->val) {
            return nullptr;
        }
        return root;
    }
}

void pruneTree_test() {
    vector<int> nums;
    TreeNode::TreeNode *root, *ans_tree;
    vector<vector<string>> ans;
    nums = {1, -1, 0, 0, 1};
    root = create_treenode(nums, true);
    ans_tree = pruneTree::pruneTree(root);
    ans = printTree::printTree(ans_tree);
    for (auto s : ans) {
        print_vector(s);
    }
    cout << "-------------" << endl;
    nums = {1, 0, 1, 0, 0, 0, 1};
    root = create_treenode(nums, true);
    ans_tree = pruneTree::pruneTree(root);
    ans = printTree::printTree(ans_tree);
    for (auto s : ans) {
        print_vector(s);
    }
    cout << "-------------" << endl;
    nums = {1, 1, 0, 1, 1, 0, 1, 0};
    root = create_treenode(nums, true);
    ans_tree = pruneTree::pruneTree(root);
    ans = printTree::printTree(ans_tree);
    for (auto s : ans) {
        print_vector(s);
    }
    cout << "-------------" << endl;
}

namespace numBusesToDestination {
    int numBusesToDestination(vector<vector<int>> &routes, int source, int target) {
        if (source == target) {
            return 0;
        }
        int n = routes.size();
        vector<vector<int>> edge(n, vector<int>(n));
        unordered_map<int, vector<int>> rec;
        for (int i = 0; i < n; i++) {
            for (int site : routes[i]) {
                for (int j : rec[site]) {
                    edge[i][j] = edge[j][i] = true;
                }
                rec[site].push_back(i);
            }
        }
        vector<int> dis(n, -1);
        queue<int> que;
        for (int bus: rec[source]) {
            dis[bus] = 1;
            que.push(bus);
        }
        while (!que.empty()) {
            int x = que.front();
            que.pop();
            for (int y = 0; y < n; ++y) {
                if (edge[x][y] && dis[y] == -1) {
                    dis[y] = dis[x] + 1;
                    que.push(y);
                }
            }
        }
        int ret = INT_MAX;
        for (int bus : rec[target]) {
            if (dis[bus] != -1) {
                ret = min(ret, dis[bus]);
            }
        }
        return ret == INT_MAX ? -1 : ret;
    }
}

void numBusesToDestination_test() {
    vector<vector<int>> routes;
    int src, tar;
    src = 1, tar = 6;
    routes = {{1, 2, 7},
              {3, 6, 7}};
    cout << numBusesToDestination::numBusesToDestination(routes, src, tar) << endl;
    src = 15, tar = 12;
    routes = {{7,  12},
              {4,  5,  15},
              {6},
              {15, 19},
              {9,  12, 13}};
    cout << numBusesToDestination::numBusesToDestination(routes, src, tar) << endl;
}

namespace ambiguousCoordinates {
    vector<string> getPos(string s) {
        vector<string> pos;
        if (s[0] != '0' || s == "0") pos.push_back(s);
        for (int p = 1; p < s.size(); ++p) {
            if ((p != 1 && s[0] == '0') || s.back() == '0') continue;
            pos.push_back(s.substr(0, p) + "." + s.substr(p));
        }
        return pos;
    }

    vector<string> ambiguousCoordinates(string s) {
        int n = s.size() - 2;
        vector<string> res;
        s = s.substr(1, s.size() - 2);
        for (int l = 1; l < n; ++l) {
            vector<string> lt = getPos(s.substr(0, l));
            if (lt.empty()) continue;
            vector<string> rt = getPos(s.substr(l));
            if (rt.empty()) continue;
            for (auto &i : lt) {
                for (auto &j : rt) {
                    res.push_back("(" + i + ", " + j + ")");
                }
            }
        }
        return res;
    }
}

void ambiguousCoordinates_test() {
    string s;
    vector<string> ans;
    s = "(123)";
    ans = ambiguousCoordinates::ambiguousCoordinates(s);
    print_vector(ans);
    s = "(00011)";
    ans = ambiguousCoordinates::ambiguousCoordinates(s);
    print_vector(ans);
    s = "(0123)";
    ans = ambiguousCoordinates::ambiguousCoordinates(s);
    print_vector(ans);
    s = "(100)";
    ans = ambiguousCoordinates::ambiguousCoordinates(s);
    print_vector(ans);
}

namespace numComponents {
    int numComponents(ListNode *head, vector<int> &nums) {
        unordered_set<int> num_set;
        for (int num : nums) {
            num_set.emplace(num);
        }
        bool in_set = false;
        int res = 0;
        while (head != nullptr) {
            if (num_set.count(head->val)) {
                if (!in_set) {
                    in_set = true;
                    res++;
                }
            } else {
                in_set = false;
            }
            head = head->next;
        }
        return res;
    }
}

void numComponents_test() {
    vector<int> lists, nums;
    ListNode *head;
    nums = {0, 1, 3};
    lists = {0, 1, 2, 3};
    head = create_nodelist(lists);
    cout << numComponents::numComponents(head, nums) << endl;
    nums = {0, 3, 1, 4};
    lists = {0, 1, 2, 3, 4};
    head = create_nodelist(lists);
    cout << numComponents::numComponents(head, nums) << endl;
    nums = {0, 1, 3, 4, 5};
    lists = {0, 1, 2, 3, 4, 5, 6};
    head = create_nodelist(lists);
    cout << numComponents::numComponents(head, nums) << endl;
}

namespace racecar {
    const int MAXN = 10001;
    const int MAXS = 14;
    int cnt[MAXN][MAXS];

    int dfs(int s, int dist) {
        if (dist == 0) {
            return 0;
        }
        if (cnt[dist][s] > 0) {
            return cnt[dist][s];
        }
        int speed = (1 << s);
        if (speed <= dist) {
            return 1 + dfs(s + 1, dist - speed);
        }
        int ans = dfs(0, speed - dist) + 2;
        for (int i = 0; i < s; ++i) {
            ans = min(ans, 2 + 2 * i + dfs(i, dist));
        }
        return cnt[dist][s] = ans;
    }

    int racecar(int target) {
        return dfs(0, target);
    }
}

void racecar_test() {
    int target;
    target = 3;
    cout << racecar::racecar(target) << endl;
    target = 6;
    cout << racecar::racecar(target) << endl;
    target = 12;
    cout << racecar::racecar(target) << endl;
};

#include<regex>

namespace mostCommonWord {
    string mostCommonWord(string paragraph, vector<string> &banned) {
        unordered_set<string> bannedSet;
        for (auto &word : banned) {
            bannedSet.emplace(word);
        }
        int maxFrequency = 0;
        unordered_map<string, int> frequencies;
        string word;
        int length = paragraph.size();
        for (int i = 0; i <= length; i++) {
            if (i < length && isalpha(paragraph[i])) {
                word.push_back(tolower(paragraph[i]));
            } else if (word.size() > 0) {
                if (!bannedSet.count(word)) {
                    frequencies[word]++;
                    maxFrequency = max(maxFrequency, frequencies[word]);
                }
                word = "";
            }
        }
        string mostCommon = "";
        for (auto &[word, frequency] : frequencies) {
            if (frequency == maxFrequency) {
                mostCommon = word;
                break;
            }
        }
        return mostCommon;
    }
}

void mostCommonWord_test() {
    string paragraph;
    vector<string> banned;
    paragraph = "Bob hit a ball, the hit BALL flew far after it was hit.";
    banned = {"hit"};
    cout << mostCommonWord::mostCommonWord(paragraph, banned) << endl;
    paragraph = "a.";
    banned = {};
    cout << mostCommonWord::mostCommonWord(paragraph, banned) << endl;
    paragraph = "W. y; v? R! n? w, U! T; V? z, M! O. z, h; s' t' t? L? G; P, o? I; k? q. Y! x, t; f! m! Y. H' W. a, Q. d! w; s; r. X; H' S? t? V; X! Z, k; j; R? v, v! H, p; Z; m! v! M; S. D; Q? P, Z. w! t? m' Y' R. c? U' z! r' T; q. J; Z? z! n? X? M; T! N, K! s? v. T? Y. e; S. q; u! V? j? O; P. U, L. m, w. w, U? o, u, z? P? P! p, I. G. x; j, I; u. q. W? S. O' N, O; F! z? s, e' W, w! y. Y' t; c; n! Z! Z; n; u! m' f, z! U! V! n' t; p? V! W; y' R; z! k, y? U; x, U; w? J; x? T; h. D. R; k; m. X' n; w' s; F! x. e! x; x; n; z. R, V, r, N' P; t; w! r? O. P. W; z; W? p. v; k; r; U; v; P. S! V; U? b' g? x, q; r! P? J' o? o. j' v; y, q, I' Y, M' i; b? l; q? W' W' f, S? s! K; G? O! m? Q. X' h. S' t! Y? Z, n. Y? V; z. r. U! R? Y! o; u, L! Z, w, d; x. z; R, o. q! y; x! l? l. S' T; R' t' y? u! V' Y; w, Z? U! k. Q? x! x; L! I! T, L? V. q? k. Z, E. K' x, O; n! Y, r. x! e! M; S; J' p' Y! N. T? i, M. u? c? x! u? P. y, P? f? S. T; V; l, y. J. s; T? Q! q; l; u; W, Z. V. h? X. w. V, n. O' s; u; j? Z' Q? g! x, I, i, r' y?";
    banned = {"y", "t", "j", "a", "s", "o", "b", "c", "d", "x", "q", "z", "l", "k", "w", "v", "r", "f", "u", "g"};
    cout << mostCommonWord::mostCommonWord(paragraph, banned) << endl;
}

namespace minimumLengthEncoding {
    int minimumLengthEncoding(vector<string> &words) {
        unordered_set<string> good(words.begin(), words.end());
        for (const string &word: words) {
            for (int k = 1; k < word.size(); ++k) {
                good.erase(word.substr(k));
            }
        }

        int ans = 0;
        for (const string &word: good) {
            ans += word.size() + 1;
        }
        return ans;
    }
}

void minimumLengthEncoding_test() {
    vector<string> words;
    words = {"time", "me", "bell"};
    cout << minimumLengthEncoding::minimumLengthEncoding(words) << endl;
    words = {"t"};
    cout << minimumLengthEncoding::minimumLengthEncoding(words) << endl;
}

namespace shortestToChar {
    vector<int> shortestToChar(string s, char c) {
        int n = s.length();
        vector<int> ans(n);
        for (int i = 0, idx = -n; i < n; ++i) {
            if (s[i] == c) {
                idx = i;
            }
            ans[i] = i - idx;
        }
        for (int i = n - 1, idx = 2 * n; i >= 0; --i) {
            if (s[i] == c) {
                idx = i;
            }
            ans[i] = min(ans[i], idx - i);
        }
        return ans;
    }
}

void shortestToChar_test() {
    string s;
    char c;
    vector<int> ans;
    s = "loveleetcode";
    c = 'e';
    ans = shortestToChar::shortestToChar(s, c);
    print_vector(ans);
    s = "aaab";
    c = 'b';
    ans = shortestToChar::shortestToChar(s, c);
    print_vector(ans);
}

namespace flipgame {
    int flipgame(vector<int> &fronts, vector<int> &backs) {
        int res = 3000, n = fronts.size();
        unordered_set<int> same;
        for (int i = 0; i < n; ++i) {
            if (fronts[i] == backs[i]) {
                same.insert(fronts[i]);
            }
        }
        for (int &x : fronts) {
            if (x < res && same.count(x) == 0) {
                res = x;
            }
        }
        for (int &x : backs) {
            if (x < res && same.count(x) == 0) {
                res = x;
            }
        }
        return res % 3000;
    }
}

void flipgame_test() {
    vector<int> fronts, backs;
    fronts = {1, 2, 4, 4, 7};
    backs = {1, 3, 4, 1, 3};
    cout << flipgame::flipgame(fronts, backs) << endl;
    fronts = {1};
    backs = {1};
    cout << flipgame::flipgame(fronts, backs) << endl;
}

namespace numFriendRequests {
    int numFriendRequests(vector<int> &ages) {
        int n = ages.size();
        sort(ages.begin(), ages.end());
        int left = 0, right = 0, ans = 0;
        for (int age : ages) {
            if (age < 15) {
                continue;
            }
            while (ages[left] <= 0.5 * age + 7) {
                ++left;
            }
            while (right + 1 < n && ages[right + 1] <= age) {
                ++right;
            }
            ans += right - left;
        }
        return ans;
    }
}

void numFriendRequests_test() {
    vector<int> ages;
    ages = {16, 16};
    cout << numFriendRequests::numFriendRequests(ages) << endl;
    ages = {16, 17, 18};
    cout << numFriendRequests::numFriendRequests(ages) << endl;
    ages = {20, 30, 100, 110, 120};
    cout << numFriendRequests::numFriendRequests(ages) << endl;
}

namespace maxProfitAssignment {
    int maxProfitAssignment(vector<int> &difficulty, vector<int> &profit, vector<int> &worker) {
        vector<pair<int, int>> jobs;
        int n = profit.size(), res = 0, i = 0, best = 0;
        for (int j = 0; j < n; ++j) {
            jobs.emplace_back(difficulty[j], profit[j]);
        }
        sort(jobs.begin(), jobs.end());
        sort(worker.begin(), worker.end());
        for (int w : worker) {
            while (i < n && w >= jobs[i].first) {
                best = max(best, jobs[i].second);
                i++;
            }
            res += best;
        }
        return res;
    }
}

void maxProfitAssignment_test() {
    vector<int> difficulty, profit, worker;
    difficulty = {2, 4, 6, 8, 10};
    profit = {10, 20, 30, 40, 50};
    worker = {4, 5, 6, 7};
    cout << maxProfitAssignment::maxProfitAssignment(difficulty, profit, worker) << endl;
    difficulty = {85, 47, 57};
    profit = {24, 66, 99};
    worker = {40, 25, 25};
    cout << maxProfitAssignment::maxProfitAssignment(difficulty, profit, worker) << endl;
}

namespace largestIsland {
    const vector<int> d = {0, -1, 0, 1, 0};

    bool valid(int n, int x, int y) {
        return x >= 0 && x < n && y >= 0 && y < n;
    }

    int dfs(const vector<vector<int>> &grid, int x, int y, vector<vector<int>> &tag, int t) {
        int n = grid.size(), res = 1;
        tag[x][y] = t;
        for (int i = 0; i < 4; i++) {
            int x1 = x + d[i], y1 = y + d[i + 1];
            if (valid(n, x1, y1) && grid[x1][y1] == 1 && tag[x1][y1] == 0) {
                res += dfs(grid, x1, y1, tag, t);
            }
        }
        return res;
    }

    int largestIsland(vector<vector<int>> &grid) {
        int n = grid.size(), res = 0;
        vector<vector<int>> tag(n, vector<int>(n));
        unordered_map<int, int> area;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (grid[i][j] == 1 && tag[i][j] == 0) {
                    int t = i * n + j + 1;
                    area[t] = dfs(grid, i, j, tag, t);
                    res = max(res, area[t]);
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (grid[i][j] == 0) {
                    int z = 1;
                    unordered_set<int> connected;
                    for (int k = 0; k < 4; k++) {
                        int x = i + d[k], y = j + d[k + 1];
                        if (!valid(n, x, y) || tag[x][y] == 0 || connected.count(tag[x][y]) > 0) {
                            continue;
                        }
                        z += area[tag[x][y]];
                        connected.insert(tag[x][y]);
                    }
                    res = max(res, z);
                }
            }
        }
        return res;
    }
}

void largestIsland_test() {
    vector<vector<int>> grid;
    grid = {{1, 0},
            {0, 1}};
    cout << largestIsland::largestIsland(grid) << endl;
    grid = {{1, 1},
            {1, 1}};
    cout << largestIsland::largestIsland(grid) << endl;
    grid = {{1, 1},
            {1, 1}};
    cout << largestIsland::largestIsland(grid) << endl;
}

namespace uniqueLetterString {
    int uniqueLetterString(string s) {
        unordered_map<char, vector<int>> index;
        for (int i = 0; i < s.size(); ++i) {
            index[s[i]].emplace_back(i);
        }
        int res = 0;
        for (auto &&[_, arr]: index) {
            arr.insert(arr.begin(), -1);
            arr.emplace_back(s.size());
            for (int i = 1; i < arr.size() - 1; i++) {
                res += (arr[i] - arr[i - 1]) * (arr[i + 1] - arr[i]);
            }
        }
        return res;
    }
}

void uniqueLetterString_test() {
    string s;
    s = "ABC";
    cout << uniqueLetterString::uniqueLetterString(s) << endl;
    s = "ABA";
    cout << uniqueLetterString::uniqueLetterString(s) << endl;
    s = "LEETCODE";
    cout << uniqueLetterString::uniqueLetterString(s) << endl;
}

namespace consecutiveNumbersSum {
    bool isKConsecutive(int n, int k) {
        if (k % 2 == 1) {
            return n % k == 0;
        } else {
            return n % k != 0 && 2 * n % k == 0;
        }
    }

    int consecutiveNumbersSum(int n) {
        int ans = 0;
        int bound = 2 * n;
        for (int k = 1; k * (k + 1) <= bound; ++k) {
            if (isKConsecutive(n, k)) {
                ans++;
            }
        }
        return ans;
    }
}

void consecutiveNumbersSum_test() {
    int n;
    n = 5;
    cout << consecutiveNumbersSum::consecutiveNumbersSum(n) << endl;
    n = 9;
    cout << consecutiveNumbersSum::consecutiveNumbersSum(n) << endl;
}

namespace largeGroupPositions {
    vector<vector<int>> largeGroupPositions(string s) {
        vector<vector<int>> ans;
        int n = s.size();
        int num = 1;
        for (int i = 0; i < n; ++i) {
            if (i == n || s[i] != s[i + 1]) {
                if (num >= 3) {
                    ans.push_back({i - num + 1, i});
                }
                num = 1;
            } else {
                num++;
            }
        }
        return ans;
    }
}

void largeGroupPositions_test() {
    string s;
    vector<vector<int>> ans;
    s = "abbxxxxzzy";
    ans = largeGroupPositions::largeGroupPositions(s);
    for (auto list : ans)
        print_vector(list);
    cout << "++++++++++" << endl;
    s = "abc";
    ans = largeGroupPositions::largeGroupPositions(s);
    for (auto list : ans)
        print_vector(list);
    cout << "++++++++++" << endl;
    s = "abcdddeeeeaabbbcd";
    ans = largeGroupPositions::largeGroupPositions(s);
    for (auto list : ans)
        print_vector(list);
    cout << "++++++++++" << endl;
    s = "aba";
    ans = largeGroupPositions::largeGroupPositions(s);
    for (auto list : ans)
        print_vector(list);
    cout << "++++++++++" << endl;
}

namespace maskPII {
    vector<string> country = {"", "+*-", "+**-", "+***-"};

    string maskPII(string s) {
        string res;
        int at = s.find("@");
        if (at != string::npos) {
            transform(s.begin(), s.end(), s.begin(), ::tolower);
            return s.substr(0, 1) + "*****" + s.substr(at - 1);
        }
        s = regex_replace(s, regex("[^0-9]"), "");
        return country[s.size() - 10] + "***-***-" + s.substr(s.size() - 4);
    }
}

void maskPII_test() {
    string s;
    s = "LeetCode@LeetCode.com";
    cout << maskPII::maskPII(s) << endl;
    s = "AB@qq.com";
    cout << maskPII::maskPII(s) << endl;
    s = "1(234)567-890";
    cout << maskPII::maskPII(s) << endl;
}

namespace flipAndInvertImage {
    vector<vector<int>> flipAndInvertImage(vector<vector<int>> &image) {
        int n = image.size();
        for (int i = 0; i < n; i++) {
            int left = 0, right = n - 1;
            while (left < right) {
                if (image[i][left] == image[i][right]) {
                    image[i][left] ^= 1;
                    image[i][right] ^= 1;
                }
                left++;
                right--;
            }
            if (left == right) {
                image[i][left] ^= 1;
            }
        }
        return image;
    }
}

void flipAndInvertImage_test() {
    vector<vector<int>> image, ans;
    image = {{1, 1, 0},
             {1, 0, 1},
             {0, 0, 0}};
    ans = flipAndInvertImage::flipAndInvertImage(image);
    for (auto row : ans) {
        print_vector(row);
    }
    cout << "+++++++++++++" << endl;
    image = {{1, 1, 0, 0},
             {1, 0, 0, 1},
             {0, 1, 1, 1},
             {1, 0, 1, 0}};
    ans = flipAndInvertImage::flipAndInvertImage(image);
    for (auto row : ans) {
        print_vector(row);
    }
    cout << "+++++++++++++" << endl;
}

namespace findReplaceString {
    string findReplaceString(string s, vector<int> &indices, vector<string> &sources, vector<string> &targets) {
        int n = s.size();
        int m = indices.size();
        vector<int> ops(m);
        iota(ops.begin(), ops.end(), 0);
        sort(ops.begin(), ops.end(), [&](int i, int j) {
            return indices[i] < indices[j];
        });
        string ans;
        int pt = 0;
        for (int i = 0; i < n;) {
            while (pt < m && indices[ops[pt]] < i) {
                ++pt;
            }
            bool succeed = false;
            while (pt < m && indices[ops[pt]] == i) {
                if (s.substr(i, sources[ops[pt]].size()) == sources[ops[pt]]) {
                    succeed = true;
                    break;
                }
                ++pt;
            }
            if (succeed) {
                ans += targets[ops[pt]];
                i += sources[ops[pt]].size();
            } else {
                ans += s[i];
                ++i;
            }
        }
        return ans;
    }
}

void findReplaceString_test() {
    string s;
    vector<int> indices;
    vector<string> sources;
    vector<string> targets;
    s = "abcd";
    indices = {0, 2};
    sources = {"a", "cd"};
    targets = {"eee", "fff"};
    cout << findReplaceString::findReplaceString(s, indices, sources, targets) << endl;
    cout << "++++++++++" << endl;
    s = "abcd";
    indices = {0, 2};
    sources = {"ab", "ec"};
    targets = {"eee", "ffff"};
    cout << findReplaceString::findReplaceString(s, indices, sources, targets) << endl;
    cout << "++++++++++" << endl;
}

namespace sumOfDistancesInTree {
    vector<int> ans, sz, dp;
    vector<vector<int>> graph;

    void dfs(int u, int f) {
        sz[u] = 1;
        dp[u] = 0;
        for (auto &v: graph[u]) {
            if (v == f) {
                continue;
            }
            dfs(v, u);
            dp[u] += dp[v] + sz[v];
            sz[u] += sz[v];
        }
    }

    void dfs2(int u, int f) {
        ans[u] = dp[u];
        for (auto &v: graph[u]) {
            if (v == f) {
                continue;
            }
            int pu = dp[u], pv = dp[v];
            int su = sz[u], sv = sz[v];

            dp[u] -= dp[v] + sz[v];
            sz[u] -= sz[v];
            dp[v] += dp[u] + sz[u];
            sz[v] += sz[u];

            dfs2(v, u);

            dp[u] = pu, dp[v] = pv;
            sz[u] = su, sz[v] = sv;
        }
    }

    vector<int> sumOfDistancesInTree(int n, vector<vector<int>> &edges) {
        ans.resize(n, 0);
        sz.resize(n, 0);
        dp.resize(n, 0);
        graph.resize(n, {});
        for (auto &edge: edges) {
            int u = edge[0], v = edge[1];
            graph[u].emplace_back(v);
            graph[v].emplace_back(u);
        }
        dfs(0, -1);
        dfs2(0, -1);
        return ans;

    }
}

void sumOfDistancesInTree_test() {
    int n = 6;
    vector<vector<int>> edges;
    vector<int> ans;
    edges = {{0, 1},
             {0, 2},
             {2, 3},
             {2, 4},
             {2, 5}};
    ans = sumOfDistancesInTree::sumOfDistancesInTree(n, edges);
    print_vector(ans);
}

namespace numSimilarGroups {
    vector<int> f;

    int find(int &x) {
        return f[x] == x ? x : f[x] = find(f[x]);
    }

    bool check(string &s1, string &s2, int length) {
        int num = 0;
        for (int i = 0; i < length; ++i) {
            if (s1[i] != s2[i]) {
                num++;
            }
            if (num > 2) {
                return false;
            }
        }
        return true;
    }

    int numSimilarGroups(vector<string> &strs) {
        int n = strs.size();
        int m = strs[0].length();
        f.resize(n);
        for (int i = 0; i < n; ++i) {
            f[i] = i;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                int fi = find(i), fj = find(j);
                if (fi == fj) {
                    continue;
                }
                if (check(strs[i], strs[j], m)) {
                    f[fi] = fj;
                }
            }
        }
        int ret = 0;
        for (int i = 0; i < n; ++i) {
            if (f[i] == i) {
                ret++;
            }
        }
        return ret;
    }
}

void numSimilarGroups_test() {
    vector<string> strs;
    strs = {"tars", "rats", "arts", "star"};
    cout << numSimilarGroups::numSimilarGroups(strs) << endl;
    strs = {"omv", "ovm"};
    cout << numSimilarGroups::numSimilarGroups(strs) << endl;
}

namespace numMagicSquaresInside {
    vector<int> m = {8, 1, 6, 7, 2, 9, 4, 3, 8, 1, 6, 7, 2, 9, 4, 3};

    bool IsMagic(vector<int> &v) {
        for (int i = 0; i < 8; i += 2)
            if (m[i] == v[0])
                return v == vector<int>(m.begin() + i, m.begin() + i + 8)
                       || v == vector<int>(m.rbegin() + 7 - i, m.rbegin() + 15 - i);
        return false;//奇数元素
    }

    int numMagicSquaresInside(vector<vector<int>> &grid) {
        int di[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
        int dj[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
        int count = 0;
        for (int i = 1; i < grid.size() - 1; i++)
            for (int j = 1; j < grid[0].size() - 1; j++)
                if (grid[i][j] == 5) {
                    vector<int> around;
                    for (int k = 0; k < 8; k++)
                        around.push_back(grid[i + di[k]][j + dj[k]]);
                    count += IsMagic(around);
                }
        return count;
    }
}

void numMagicSquaresInside_test() {
    vector<vector<int>> grid;
    grid = {{4, 3, 8, 4},
            {9, 5, 1, 9},
            {2, 7, 6, 2}};
    cout << numMagicSquaresInside::numMagicSquaresInside(grid) << endl;
    grid = {{8}};
    cout << numMagicSquaresInside::numMagicSquaresInside(grid) << endl;
}

namespace canVisitAllRooms {
    vector<int> vis;
    int num;

    void dfs(vector<vector<int>> &rooms, int x) {
        vis[x] = true;
        num++;
        for (auto &it : rooms[x]) {
            if (!vis[it]) {
                dfs(rooms, it);
            }
        }
    }

    bool canVisitAllRooms(vector<vector<int>> &rooms) {
        int n = rooms.size();
        num = 0;
        vis.resize(n);
        dfs(rooms, 0);
        return num == n;
    }
}

void canVisitAllRooms_test() {
    vector<vector<int>> rooms;
    rooms = {{1},
             {2},
             {3},
             {}};
    cout << canVisitAllRooms::canVisitAllRooms(rooms) << endl;
    rooms = {{1, 3},
             {3, 0, 1},
             {2},
             {0}};
    cout << canVisitAllRooms::canVisitAllRooms(rooms) << endl;
}

namespace backspaceCompare {
    bool backspaceCompare(string s, string t) {
        stack<char> stks, stkt;
        for (auto c : s) {
            if (c == '#') {
                if (stks.empty()) {
                    continue;
                }
                stks.pop();
            } else {
                stks.push(c);
            }
        }
        for (auto c : t) {
            if (c == '#') {
                if (stkt.empty()) {
                    continue;
                }
                stkt.pop();
            } else {
                stkt.push(c);
            }
        }
        while (!stks.empty() && !stkt.empty()) {
            auto cs = stks.top();
            auto ct = stkt.top();
            stks.pop();
            stkt.pop();
            if (cs != ct) {
                return false;
            }
        }
        return stks.empty() && stkt.empty() ? true : false;
    }
}

void backspaceCompare_test() {
    string s, t;
    s = "ab#c";
    t = "ad#c";
    cout << backspaceCompare::backspaceCompare(s, t) << endl;
    s = "ab##";
    t = "c#d#";
    cout << backspaceCompare::backspaceCompare(s, t) << endl;
    s = "a#c";
    t = "b";
    cout << backspaceCompare::backspaceCompare(s, t) << endl;
}

namespace longestMountain {
    int longestMountain(vector<int> &arr) {
        int n = arr.size();
        if (n == 0) {
            return 0;
        }
        vector<int> left(n), right(n);
        for (int i = 1; i < n; ++i) {
            left[i] = (arr[i - 1] < arr[i] ? left[i - 1] + 1 : 0);
        }
        for (int i = n - 2; i >= 0; --i) {
            right[i] = (arr[i + 1] < arr[i] ? right[i + 1] + 1 : 0);
        }
        int ans = 0;
        for (int i = 0; i < n; ++i) {
            if (left[i] > 0 && right[i] > 0) {
                ans = max(ans, left[i] + right[i] + 1);
            }
        }
        return ans;
    }
}

void longestMountain_test() {
    vector<int> arr;
    arr = {2, 1, 4, 7, 3, 2, 5};
    cout << longestMountain::longestMountain(arr) << endl;
    arr = {2, 2, 2};
    cout << longestMountain::longestMountain(arr) << endl;
}

namespace isNStraightHand {
    bool isNStraightHand(vector<int> &hand, int groupSize) {
        int n = hand.size();
        if (n % groupSize != 0) {
            return false;
        }
        sort(hand.begin(), hand.end());
        unordered_map<int, int> cnt;
        for (auto &num : hand) {
            cnt[num]++;
        }
        for (auto &x : hand) {
            if (!cnt.count(x)) {
                continue;
            }
            for (int j = 0; j < groupSize; j++) {
                int num = x + j;
                if (!cnt.count(num)) {
                    return false;
                }
                cnt[num]--;
                if (cnt[num] == 0) {
                    cnt.erase(num);
                }
            }
        }
        return true;
    }
}

void isNStraightHand_test() {
    vector<int> hand;
    int group_size;
    hand = {1, 2, 3, 6, 2, 3, 4, 7, 8};
    group_size = 3;
    cout << isNStraightHand::isNStraightHand(hand, group_size) << endl;
    hand = {1, 2, 3, 4, 5};
    group_size = 4;
    cout << isNStraightHand::isNStraightHand(hand, group_size) << endl;
}

namespace shortestPathLength {
    int shortestPathLength(vector<vector<int>> &graph) {
        int n = graph.size();
        queue<tuple<int, int, int>> q;
        vector<vector<int>> seen(n, vector<int>(1 << n));
        for (int i = 0; i < n; ++i) {
            q.emplace(i, 1 << i, 0);
            seen[i][1 << i] = true;
        }

        int ans = 0;
        while (!q.empty()) {
            auto[u, mask, dist] = q.front();
            q.pop();
            if (mask == (1 << n) - 1) {
                ans = dist;
                break;
            }
            // 搜索相邻的节点
            for (int v: graph[u]) {
                // 将 mask 的第 v 位置为 1
                int mask_v = mask | (1 << v);
                if (!seen[v][mask_v]) {
                    q.emplace(v, mask_v, dist + 1);
                    seen[v][mask_v] = true;
                }
            }
        }
        return ans;
    }
}

void shortestPathLength_test() {
    vector<vector<int>> graph;
    graph = {{1, 2, 3},
             {0},
             {0},
             {0}};
    cout << shortestPathLength::shortestPathLength(graph) << endl;
    graph = {{1},
             {0, 2, 4},
             {1, 3, 4},
             {2},
             {1, 2}};
    cout << shortestPathLength::shortestPathLength(graph) << endl;
}

namespace maxDistToClosest {
    int maxDistToClosest(vector<int> &seats) {
        int l = 0;
        int res = 0;
        while (l < seats.size() && seats[l] == 0) {
            ++l;
        }
        res = max(res, l);
        while (l < seats.size()) {
            int r = l + 1;
            while (r < seats.size() && seats[r] == 0) {
                ++r;
            }
            if (r == seats.size()) {
                res = max(res, r - l - 1);
            } else {
                res = max(res, (r - l) / 2);
            }
            l = r;
        }
        return res;
    }
}

void maxDistToClosest_test() {
    vector<int> seats;
//    seats = {1, 0, 0, 0, 1, 0, 1};
//    cout << maxDistToClosest::maxDistToClosest(seats) << endl;
//    seats = {1, 0, 0, 0};
//    cout << maxDistToClosest::maxDistToClosest(seats) << endl;
//    seats = {0, 1};
//    cout << maxDistToClosest::maxDistToClosest(seats) << endl;
    seats = {0, 0, 1, 1};
    cout << maxDistToClosest::maxDistToClosest(seats) << endl;
}

namespace rectangleArea {
    int rectangleArea(vector<vector<int>> &rectangles) {
        int n = rectangles.size();
        vector<int> hbound;
        for (const auto &rect : rectangles) {
            hbound.push_back(rect[1]);
            hbound.push_back(rect[3]);
        }
        sort(hbound.begin(), hbound.end());
        hbound.erase(unique(hbound.begin(), hbound.end()), hbound.end());
        int m = hbound.size();
        vector<int> seg(m - 1);
        vector<tuple<int, int, int>> sweep;
        for (int i = 0; i < n; ++i) {
            // 左边界
            sweep.emplace_back(rectangles[i][0], i, 1);
            // 右边界
            sweep.emplace_back(rectangles[i][2], i, -1);
        }
        sort(sweep.begin(), sweep.end());

        long long ans = 0;
        for (int i = 0; i < sweep.size(); ++i) {
            int j = i;
            while (j + 1 < sweep.size() && get<0>(sweep[i]) == get<0>(sweep[j + 1])) {
                ++j;
            }
            if (j + 1 == sweep.size()) {
                break;
            }
            // 一次性地处理掉一批横坐标相同的左右边界
            for (int k = i; k <= j; ++k) {
                auto&&[_, idx, diff] = sweep[k];
                int left = rectangles[idx][1], right = rectangles[idx][3];
                for (int x = 0; x < m - 1; ++x) {
                    if (left <= hbound[x] && hbound[x + 1] <= right) {
                        seg[x] += diff;
                    }
                }
            }
            int cover = 0;
            for (int k = 0; k < m - 1; ++k) {
                if (seg[k] > 0) {
                    cover += (hbound[k + 1] - hbound[k]);
                }
            }
            ans += static_cast<long long>(cover) * (get<0>(sweep[j + 1]) - get<0>(sweep[j]));
            i = j;
        }
        return ans % static_cast<int>(1e9 + 7);
    }
}

void rectangleArea_test() {
    vector<vector<int>> rectangles;
    rectangles = {{0, 0, 2, 2},
                  {1, 0, 2, 3},
                  {1, 0, 3, 1}};
    cout << rectangleArea::rectangleArea(rectangles) << endl;
    rectangles = {{0, 0, 1000000000, 1000000000}};
    cout << rectangleArea::rectangleArea(rectangles) << endl;
}

namespace loudAndRich {
    vector<int> loudAndRich(vector<vector<int>> &richer, vector<int> &quiet) {
        int n = quiet.size();
        vector<vector<int>> g(n);
        for (auto &r : richer) {
            g[r[1]].emplace_back(r[0]);
        }
        vector<int> ans(n, -1);
        function<void(int)> dfs = [&](int x) {
            if (ans[x] != -1) {
                return;
            }
            ans[x] = x;
            for (int y : g[x]) {
                dfs(y);
                if (quiet[ans[y]] < quiet[ans[x]]) {
                    ans[x] = ans[y];
                }
            }
        };
        for (int i = 0; i < n; ++i) {
            dfs(i);
        }
        return ans;
    }
}

void loudAndRich_test() {
    vector<vector<int>> richer;
    vector<int> quiet;
    richer = {{1, 0},
              {2, 1},
              {3, 1},
              {3, 7},
              {4, 3},
              {5, 3},
              {6, 3}};
    quiet = {3, 2, 5, 4, 6, 1, 7, 0};
    vector<int> ans;
    ans = loudAndRich::loudAndRich(richer, quiet);
    print_vector(ans);
    richer = {};
    quiet = {0};
    ans = loudAndRich::loudAndRich(richer, quiet);
    print_vector(ans);
}

namespace peakIndexInMountainArray {
    int peakIndexInMountainArray(vector<int> &arr) {
        int n = arr.size();
        int left = 1, right = n - 2, ans = 0;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (arr[mid] > arr[mid + 1]) {
                ans = mid;
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        return ans;
    }
}

void peakIndexInMountainArray_test() {
    vector<int> arr;
    arr = {0, 1, 0};
    cout << peakIndexInMountainArray::peakIndexInMountainArray(arr) << endl;
    arr = {0, 2, 1, 0};
    cout << peakIndexInMountainArray::peakIndexInMountainArray(arr) << endl;
    arr = {0, 10, 5, 2};
    cout << peakIndexInMountainArray::peakIndexInMountainArray(arr) << endl;
}

namespace decodeAtIndex {
    string decodeAtIndex(string s, int k) {
        long size = 0;
        int N = s.size();

        // Find size = length of decoded string
        for (int i = 0; i < N; ++i) {
            if (isdigit(s[i]))
                size *= s[i] - '0';
            else
                size++;
        }
        for (int i = N - 1; i >= 0; --i) {
            k %= size;
            if (k == 0 && isalpha(s[i]))
                return (string) "" + s[i];

            if (isdigit(s[i]))
                size /= s[i] - '0';
            else
                size--;
        }
        return "";
    }
}

void decodeAtIndex_test() {
    string s;
    int k;
//    s = "leet2code3";
//    k = 10;
//    cout << decodeAtIndex::decodeAtIndex(s, k) << endl;
//    s = "ha22";
//    k = 5;
//    cout << decodeAtIndex::decodeAtIndex(s, k) << endl;
    s = "a2345678999999999999999";
    k = 1;
    cout << decodeAtIndex::decodeAtIndex(s, k) << endl;
}

namespace numRescueBoats {
    int numRescueBoats(vector<int> &people, int limit) {
        int ans = 0;
        sort(people.begin(), people.end());
        int light = 0, heavy = people.size() - 1;
        while (light <= heavy) {
            if (people[light] + people[heavy] > limit) {
                --heavy;
            } else {
                ++light;
                --heavy;
            }
            ++ans;
        }
        return ans;
    }
}

void numRescueBoats_test() {
    vector<int> people;
    int limit;
//    people = {1, 2};
//    limit = 3;
//    cout << numRescueBoats::numRescueBoats(people, limit) << endl;
    people = {3, 2, 2, 1};
    limit = 3;
    cout << numRescueBoats::numRescueBoats(people, limit) << endl;
    people = {3, 5, 3, 4};
    limit = 5;
    cout << numRescueBoats::numRescueBoats(people, limit) << endl;
}

namespace reachableNodes {
    int encode(int u, int v, int n) {
        return u * n + v;
    }

    int reachableNodes(vector<vector<int>> &edges, int maxMoves, int n) {
        vector<vector<pair<int, int>>> adList(n);
        for (auto &edge: edges) {
            int u = edge[0], v = edge[1], nodes = edge[2];
            adList[u].emplace_back(v, nodes);
            adList[v].emplace_back(u, nodes);
        }
        unordered_map<int, int> used;
        unordered_set<int> visited;
        int reachableNodes = 0;
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
        pq.emplace(0, 0);
        while (!pq.empty() && pq.top().first <= maxMoves) {
            auto[step, u] = pq.top();
            pq.pop();
            if (visited.count(u)) {
                continue;
            }
            visited.emplace(u);
            reachableNodes++;
            for (auto[v, nodes] : adList[u]) {
                if (nodes + step + 1 <= maxMoves && !visited.count(v)) {
                    pq.emplace(nodes + step + 1, v);
                }
                used[encode(u, v, n)] = min(nodes, maxMoves - step);
            }
        }
        for (auto &edge : edges) {
            int u = edge[0], v = edge[1], nodes = edge[2];
            reachableNodes += min(nodes, used[encode(u, v, n)] + used[encode(v, u, n)]);
        }
        return reachableNodes;
    }
}

void reachableNodes_test() {
    vector<vector<int>> edges;
    int n, max_moves;
    edges = {{0, 1, 10},
             {0, 2, 1},
             {1, 2, 2}};
    max_moves = 6;
    cout << reachableNodes::reachableNodes(edges, max_moves, n) << endl;
}

namespace uncommonFromSentences {
    vector<string> wordsSplit(string str) {
        vector<string> words;
        size_t start = 0, end;
        while ((end = str.find(' ', start)) != std::string::npos) {
            words.push_back(str.substr(start, end - start));  // 提取子字符串
            start = end + 1;  // 更新开始位置
        }
        words.push_back(str.substr(start));
        return words;
    }

    vector<string> uncommonFromSentences(string s1, string s2) {
        unordered_map<string, int> map;
        vector<string> words1, words2, ans;
        words1 = uncommonFromSentences::wordsSplit(s1);
        words2 = uncommonFromSentences::wordsSplit(s2);
        for (auto &word : words1) {
            map[word]++;
        }
        for (auto &word : words2) {
            map[word]++;
        }
        for (unordered_map<string, int>::const_iterator it = map.cbegin(); it != map.cend(); ++it) {
            if (it->second == 1) {
                ans.push_back(it->first);
            }
        }
        return ans;
    }
}

void uncommonFromSentences_test() {
    string s1, s2;
    vector<string> ans;
    s1 = "this apple is sweet";
    s2 = "this apple is sour";
    ans = uncommonFromSentences::uncommonFromSentences(s1, s2);
    print_vector(ans);
    s1 = "apple apple";
    s2 = "banana";
    ans = uncommonFromSentences::uncommonFromSentences(s1, s2);
    print_vector(ans);
}

namespace spiralMatrixIII {
    vector<vector<int>> spiralMatrixIII(int rows, int cols, int rStart, int cStart) {
        vector<vector<int>> res;
        vector<pair<int, int>> around = {{0,  1},
                                         {1,  0},
                                         {0,  -1},
                                         {-1, 0}};  //顺时针方向
        int x = rStart, y = cStart, num = 1, dir = 0;  //{x, y}为当前位置，num为当前查找的数字，dir为当前的方向
        int Left = cStart - 1, Right = cStart + 1, Upper = rStart - 1, Bottom = rStart + 1;  //四个方向的边界
        while (num <= rows * cols) {
            if (x >= 0 && x < rows && y >= 0 && y < cols) {  //{x， y}位置在矩阵中
                res.push_back({x, y});
                num += 1;
            }
            if (dir == 0 && y == Right) {  //向右到右边界
                dir += 1;  //调转方向向下
                Right += 1;  //右边界右移
            } else if (dir == 1 && x == Bottom) {  //向下到底边界
                dir += 1;
                Bottom += 1;  //底边界下移
            } else if (dir == 2 && y == Left) {  //向左到左边界
                dir += 1;
                Left--;  //左边界左移
            } else if (dir == 3 && x == Upper) {  //向上到上边界
                dir = 0;
                Upper--;  //上边界上移
            }
            x += around[dir].first;   //下一个节点
            y += around[dir].second;
        }
        return res;
    }
}

void spiralMatrixIII_test() {
    int rows, cols, rStart, cStart;
    vector<vector<int>> ans;
    rows = 1, cols = 4;
    rStart = 0, cStart = 0;
    ans = spiralMatrixIII::spiralMatrixIII(rows, cols, rStart, cStart);
    for (auto list : ans) {
        print_vector(list);
    }
    cout << "------------" << endl;
    rows = 5, cols = 6, rStart = 1, cStart = 4;
    ans = spiralMatrixIII::spiralMatrixIII(rows, cols, rStart, cStart);
    for (auto list : ans) {
        print_vector(list);
    }
    cout << "------------" << endl;
}

namespace possibleBipartition {
    bool dfs(int curnode, int nowcolor, vector<int> &color, const vector<vector<int>> &g) {
        color[curnode] = nowcolor;
        for (auto &nextnode : g[curnode]) {
            if (color[nextnode] && color[nextnode] == color[curnode]) {
                return false;
            }
            if (!color[nextnode] && !dfs(nextnode, 3 ^ nowcolor, color, g)) {
                return false;
            }
        }
        return true;
    }

    bool possibleBipartition(int n, vector<vector<int>> &dislikes) {
        vector<int> color(n + 1, 0);
        vector<vector<int>> g(n + 1);
        for (auto &p : dislikes) {
            g[p[0]].push_back(p[1]);
            g[p[1]].push_back(p[0]);
        }
        for (int i = 1; i <= n; ++i) {
            if (color[i] == 0 && !dfs(i, 1, color, g)) {
                return false;
            }
        }
        return true;
    }
}

void possibleBipartition_test() {
    vector<vector<int>> dislikes;
    int n;
    n = 4;
    dislikes = {{1, 2},
                {1, 3},
                {2, 4}};
    cout << possibleBipartition::possibleBipartition(n, dislikes) << endl;
    n = 3;
    dislikes = {{1, 2},
                {1, 3},
                {2, 3}};
    cout << possibleBipartition::possibleBipartition(n, dislikes) << endl;
    n = 5;
    dislikes = {{1, 2},
                {2, 3},
                {3, 4},
                {4, 5},
                {1, 5}};
    cout << possibleBipartition::possibleBipartition(n, dislikes) << endl;
}

namespace superEggDrop {
    unordered_map<int, int> memo;

    int dp(int k, int n) {
        if (memo.find(n * 100 + k) == memo.end()) {
            int ans;
            if (n == 0) {
                ans = 0;
            } else if (k == 1) {
                ans = n;
            } else {
                int lo = 1, hi = n;
                while (lo + 1 < hi) {
                    int x = (lo + hi) / 2;
                    int t1 = dp(k - 1, x - 1);
                    int t2 = dp(k, n - x);

                    if (t1 < t2) {
                        lo = x;
                    } else if (t1 > t2) {
                        hi = x;
                    } else {
                        lo = hi = x;
                    }
                }

                ans = 1 + min(max(dp(k - 1, lo - 1), dp(k, n - lo)),
                              max(dp(k - 1, hi - 1), dp(k, n - hi)));
            }

            memo[n * 100 + k] = ans;
        }

        return memo[n * 100 + k];
    }

    int superEggDrop(int k, int n) {
        return dp(k, n);
    }
}

void superEggDrop_test() {
    int k, n;
    k = 1;
    n = 2;
    cout << superEggDrop::superEggDrop(k, n) << endl;
    k = 2;
    n = 4;
    cout << superEggDrop::superEggDrop(k, n) << endl;
    k = 3;
    n = 14;
    cout << superEggDrop::superEggDrop(k, n) << endl;
}

namespace fairCandySwap {
    vector<int> fairCandySwap(vector<int> &aliceSizes, vector<int> &bobSizes) {
        unordered_map<int, int> map_alice, map_bob;
        int alice_sum = accumulate(aliceSizes.begin(), aliceSizes.end(), 0);
        int bob_sum = accumulate(bobSizes.begin(), bobSizes.end(), 0);
        unordered_set<int> alice_set, bob_set;
        for (auto &num : aliceSizes) {
            map_alice[num] = (bob_sum - alice_sum + 2 * num) / 2;
        }
        for (auto &num : bobSizes) {
            map_bob[num] = (alice_sum - bob_sum + 2 * num) / 2;
        }
        for (auto &num : aliceSizes) {
            if (map_bob.find(map_alice[num]) != map_bob.end()) {
                return {num, map_alice[num]};
            }
        }
    }
}

void fairCandySwap_test() {
    vector<int> aliceSizes, bobSizes, ans;
    aliceSizes = {1, 1}, bobSizes = {2, 2};
    ans = fairCandySwap::fairCandySwap(aliceSizes, bobSizes);
    print_vector(ans);
    aliceSizes = {1, 2}, bobSizes = {2, 3};
    ans = fairCandySwap::fairCandySwap(aliceSizes, bobSizes);
    print_vector(ans);
    aliceSizes = {2}, bobSizes = {1, 3};
    ans = fairCandySwap::fairCandySwap(aliceSizes, bobSizes);
    print_vector(ans);
    aliceSizes = {1, 2, 5}, bobSizes = {2, 4};
    ans = fairCandySwap::fairCandySwap(aliceSizes, bobSizes);
    print_vector(ans);
}

namespace constructFromPrePost {
    TreeNode::TreeNode *constructFromPrePost(vector<int> &preorder, vector<int> &postorder) {
        int n = preorder.size();
        unordered_map<int, int> postMap;
        for (int i = 0; i < n; i++) {
            postMap[postorder[i]] = i;
        }
        function<TreeNode::TreeNode *(int, int, int, int)> dfs = [&](int preLeft, int preRight, int postLeft,
                                                                     int postRight) -> TreeNode::TreeNode * {
            if (preLeft > preRight) {
                return nullptr;
            }
            int leftCount = 0;
            if (preLeft < preRight) {
                leftCount = postMap[preorder[preLeft + 1]] - postLeft + 1;
            }
            return new TreeNode::TreeNode(preorder[preLeft],
                                          dfs(preLeft + 1, preLeft + leftCount, postLeft, postLeft + leftCount - 1),
                                          dfs(preLeft + leftCount + 1, preRight, postLeft + leftCount, postRight - 1));
        };
        return dfs(0, n - 1, 0, n - 1);
    }
}

void constructFromPrePost_test() {
    vector<int> preorder, postorder;
    TreeNode::TreeNode *ans;
    preorder = {1, 2, 4, 5, 3, 6, 7}, postorder = {4, 5, 2, 6, 7, 3, 1};
    ans = constructFromPrePost::constructFromPrePost(preorder, postorder);
    printTree::printTree(ans);
    preorder = {1}, postorder = {1};
    ans = constructFromPrePost::constructFromPrePost(preorder, postorder);
    printTree::printTree(ans);
}

namespace findAndReplacePattern {
    bool match(string &word, string &pattern) {
        unordered_map<char, char> map;
        for (int i = 0; i < word.size(); ++i) {
            char x = word[i], y = pattern[i];
            if (!map.count(x)) {
                map[x] = y;
            } else if (map[x] != y) {
                return false;
            }
        }
        return true;
    }

    vector<string> findAndReplacePattern(vector<string> &words, string pattern) {
        vector<string> ans;
        for (auto &word : words) {
            if (match(word, pattern) && match(pattern, word)) {
                ans.emplace_back(word);
            }
        }
        return ans;
    }
}

void findAndReplacePattern_test() {
    vector<string> words, ans;
    string pattern;
    words = {"abc", "deq", "mee", "aqq", "dkd", "ccc"};
    pattern = "abb";
    ans = findAndReplacePattern::findAndReplacePattern(words, pattern);
    print_vector(ans);
    words = {"a", "b", "c"};
    pattern = "a";
    ans = findAndReplacePattern::findAndReplacePattern(words, pattern);
    print_vector(ans);
}

namespace sumSubseqWidths {
    int sumSubseqWidths(vector<int> &nums) {
        sort(nums.begin(), nums.end());
        long long res = 0, mod = 1e9 + 7;
        long long x = nums[0], y = 2;
        for (int j = 1; j < nums.size(); ++j) {
            res = (res + nums[j] * (y - 1) - x) % mod;
            x = (x * 2 + nums[j]) % mod;
            y = y * 2 % mod;
        }
        return (res + mod) % mod;
    }
}

void sumSubseqWidths_test() {
    vector<int> nums;
    nums = {2, 1, 3};
    cout << sumSubseqWidths::sumSubseqWidths(nums) << endl;
    nums = {2};
    cout << sumSubseqWidths::sumSubseqWidths(nums) << endl;
}

namespace surfaceArea {
    int surfaceArea(vector<vector<int>> &grid) {
        int dr[]{0, 1, 0, -1};
        int dc[]{1, 0, -1, 0};

        int N = grid.size();
        int ans = 0;

        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) {
                if (grid[r][c] > 0) {
                    ans += 2;
                    for (int k = 0; k < 4; ++k) {
                        int nr = r + dr[k];
                        int nc = c + dc[k];
                        int nv = 0;
                        if (0 <= nr && nr < N && 0 <= nc && nc < N) {
                            nv = grid[nr][nc];
                        }

                        ans += max(grid[r][c] - nv, 0);
                    }
                }
            }
        }
        return ans;
    }
}

void surfaceArea_test() {
    vector<vector<int>> grid;
    grid = {{1, 2},
            {3, 4}};
    cout << surfaceArea::surfaceArea(grid) << endl;
    grid = {{1, 1, 1},
            {1, 0, 1},
            {1, 1, 1}};
    cout << surfaceArea::surfaceArea(grid) << endl;
}

namespace numSpecialEquivGroups {
    int numSpecialEquivGroups(vector<string> &words) {
        unordered_set<string> set;
        for (auto &s : words) {
            char count[52];
            for (int i = 0; i < s.size(); ++i) {
                count[s[i] - 'a' + 26 * (i % 2)]++;
            }
            set.insert(string(count));
        }
        return set.size();
    }
}

void numSpecialEquivGroups_test() {
    vector<string> words;
    words = {"abcd", "cdab", "cbad", "xyzz", "zzxy", "zzyx"};
    cout << numSpecialEquivGroups::numSpecialEquivGroups(words) << endl;
    words = {"abc", "acb", "bac", "bca", "cab", "cba"};
    cout << numSpecialEquivGroups::numSpecialEquivGroups(words) << endl;
}

namespace allPossibleFBT {
    vector<TreeNode::TreeNode *> allPossibleFBT(int n) {
        vector<TreeNode::TreeNode *> fullBinaryTrees;
        if (n % 2 == 0) {
            return fullBinaryTrees;
        }
        if (n == 1) {
            fullBinaryTrees = {new TreeNode::TreeNode(1)};
            return fullBinaryTrees;
        }
        for (int i = 1; i < n; i += 2) {
            vector<TreeNode::TreeNode *> leftSubtrees = allPossibleFBT(i);
            vector<TreeNode::TreeNode *> rightSubtrees = allPossibleFBT(n - i - 1);
            for (TreeNode::TreeNode *leftSubtree : leftSubtrees) {
                for (TreeNode::TreeNode *rightSubtree : rightSubtrees) {
                    TreeNode::TreeNode *root = new TreeNode::TreeNode(1, leftSubtree, rightSubtree);
                    fullBinaryTrees.emplace_back(root);
                }
            }
        }
        return fullBinaryTrees;
    }
}

void allPossibleFBT_test() {
    int n;
    vector<TreeNode::TreeNode *> ans;
    n = 7;
    ans = allPossibleFBT::allPossibleFBT(n);
    for (auto *node : ans) {
        auto s = TreeNode::print_tree(node);
        cout << s << endl;
    }
    cout << "-----------" << endl;
    n = 3;
    ans = allPossibleFBT::allPossibleFBT(n);
    for (auto *node : ans) {
        auto s = TreeNode::print_tree(node);
        cout << s << endl;
    }
}

namespace isMonotonic {
    bool isMonotonic(vector<int> &nums) {
        int mono = -1;
        if (nums.size() == 1) {
            return true;
        }
        for (int i = 1; i < nums.size(); ++i) {
            if (nums[i - 1] < nums[i]) {
                if (mono == -1) {
                    mono = 1;
                } else if (mono == 2) {
                    return false;
                }
            } else if (nums[i - 1] > nums[i]) {
                if (mono == -1) {
                    mono = 2;
                } else if (mono == 1) {
                    return false;
                }
            }
        }
        return true;
    }
}

void isMonotonic_test() {
    vector<int> nums;
    nums = {1, 2, 2,};
    cout << isMonotonic::isMonotonic(nums) << endl;
    nums = {6, 5, 4, 4};
    cout << isMonotonic::isMonotonic(nums) << endl;
    nums = {1, 3, 2};
    cout << isMonotonic::isMonotonic(nums) << endl;
}

namespace increasingBST {
    TreeNode::TreeNode *resNode;

    void inorder(TreeNode::TreeNode *node) {
        if (node == nullptr) {
            return;
        }
        inorder(node->left);

        // 在中序遍历的过程中修改节点指向
        resNode->right = node;
        node->left = nullptr;
        resNode = node;

        inorder(node->right);
    }

    TreeNode::TreeNode *increasingBST(TreeNode::TreeNode *root) {
        TreeNode::TreeNode *dummyNode = new TreeNode::TreeNode(-1);
        resNode = dummyNode;
        inorder(root);
        return dummyNode->right;
    }
}

void increasingBST_test() {
    vector<int> nums;
    TreeNode::TreeNode *root, *ans;
    vector<vector<string>> output;
    string s;
    nums = {5, 3, 6, 2, 4, 0, 8, 1, 0, 0, 0, 7, 9};
    root = create_treenode(nums, false);
    ans = increasingBST::increasingBST(root);
    s = TreeNode::print_tree(ans);
    cout << s << endl;
    cout << "-----------" << endl;
    nums = {5, 1, 7};
    root = create_treenode(nums, false);
    ans = increasingBST::increasingBST(root);
    s = TreeNode::print_tree(ans);
    cout << s << endl;
    cout << "-----------" << endl;
}

namespace subarrayBitwiseORs {
    int subarrayBitwiseORs(vector<int> &arr) {
        // ors 保留前面子数组的或运算的所有结果值
        unordered_set<int> res, ors;
        for (int x:arr) {
            unordered_set<int> tmp;
            // 将 ors 中的各元素与当前元素进行或运算
            for (auto it = ors.begin(); it != ors.end(); it++)
                tmp.insert((*it) | x);
            // 插入当前元素
            tmp.insert(x);
            ors = tmp;
            // 与原来保存所有或运算结果值的res做并集
            for (auto it = tmp.begin(); it != tmp.end(); it++)
                res.insert(*it);
        }
        return res.size();
    }
}

void subarrayBitwiseORs_test() {
    vector<int> arr;
    arr = {0};
    cout << subarrayBitwiseORs::subarrayBitwiseORs(arr) << endl;
    arr = {1, 1, 2};
    cout << subarrayBitwiseORs::subarrayBitwiseORs(arr) << endl;
    arr = {1, 2, 4};
    cout << subarrayBitwiseORs::subarrayBitwiseORs(arr) << endl;
}

namespace orderlyQueue {
    string orderlyQueue(string s, int k) {
        if (k == 1) {
            string smallest = s;
            int n = s.size();
            for (int i = 1; i < n; i++) {
                char c = s[0];
                s = s.substr(1);
                s.push_back(c);
                if (s < smallest) {
                    smallest = s;
                }
            }
            return smallest;
        } else {
            sort(s.begin(), s.end());
            return s;
        }
    }
}

void orderlyQueue_test() {
    string s;
    int k;
    s = "cba";
    k = 1;
    cout << orderlyQueue::orderlyQueue(s, k) << endl;
    s = "baaca";
    k = 3;
    cout << orderlyQueue::orderlyQueue(s, k) << endl;
}

namespace atMostNGivenDigitSet {
    int atMostNGivenDigitSet(vector<string> &digits, int n) {
        string s = to_string(n);
        int m = digits.size(), k = s.size();
        vector<vector<int>> dp(k + 1, vector<int>(2));
        dp[0][1] = 1;
        for (int i = 1; i <= k; ++i) {
            for (int j = 0; j < m; ++j) {
                if (digits[j][0] == s[i - 1]) {
                    dp[i][1] = dp[i - 1][1];
                } else if (digits[j][0] < s[i - 1]) {
                    dp[i][0] += dp[i - 1][1];
                } else {
                    break;
                }
            }
            if (i > 1) {
                dp[i][0] += m + dp[i - 1][0] * m;
            }
        }
        return dp[k][0] + dp[k][1];
    }
}

void atMostNGivenDigitSet_test() {
    vector<string> digits;
    int n;
    digits = {"1", "3", "5", "7"};
    n = 100;
    cout << atMostNGivenDigitSet::atMostNGivenDigitSet(digits, n) << endl;
    digits = {"1", "4", "9"};
    n = 1000000000;
    cout << atMostNGivenDigitSet::atMostNGivenDigitSet(digits, n) << endl;
    digits = {"7"};
    n = 8;
    cout << atMostNGivenDigitSet::atMostNGivenDigitSet(digits, n) << endl;
}

namespace numPermsDISequence {
    int numPermsDISequence(string s) {
        int i, j, size = s.size(), sum = 0, mod = 1000000007;
        vector<vector<int>> dp(size + 1, vector<int>(size + 1));
        dp[0][0] = 1;
        for (i = 1; i <= size; ++i) {
            if (s[i - 1] == 'D') {
                dp[i][i] = 0;
                for (j = i - 1; j >= 0; --j) {
                    dp[i][j] = (dp[i][j + 1] + dp[i - 1][j]) % mod;
                }
            } else {
                dp[i][0] = 0;
                for (j = 1; j <= i; ++j) {
                    dp[i][j] = (dp[i][j - 1] + dp[i - 1][j - 1]) % mod;
                }
            }
        }
        for (j = 0; j <= size; ++j) {
            sum = (sum + dp[size][j]) % mod;
        }
        return sum;
    }
}

void numPermsDISequence_test() {
    string s;
    s = "DID";
    cout << numPermsDISequence::numPermsDISequence(s) << endl;
    s = "D";
    cout << numPermsDISequence::numPermsDISequence(s) << endl;
}

namespace totalFruit {
    int totalFruit(vector<int> &fruits) {
        int n = fruits.size();
        unordered_map<int, int> cnt;

        int left = 0, ans = 0;
        for (int right = 0; right < n; ++right) {
            ++cnt[fruits[right]];
            while (cnt.size() > 2) {
                auto it = cnt.find(fruits[left]);
                --it->second;
                if (it->second == 0) {
                    cnt.erase(it);
                }
                ++left;
            }
            ans = max(ans, right - left + 1);
        }
        return ans;
    }
}

void totalFruit_test() {
    vector<int> fruits;
    fruits = {1, 2, 1};
    cout << totalFruit::totalFruit(fruits) << endl;
    fruits = {0, 1, 2, 2};
    cout << totalFruit::totalFruit(fruits) << endl;
    fruits = {1, 2, 3, 2, 2};
    cout << totalFruit::totalFruit(fruits) << endl;
    fruits = {3, 3, 3, 1, 2, 1, 1, 2, 3, 3, 4};
    cout << totalFruit::totalFruit(fruits) << endl;
}

namespace sortArrayByParity {
    vector<int> sortArrayByParity(vector<int> &nums) {
        int left = 0, right = nums.size() - 1;
        while (left < right) {
            while (left < right and nums[left] % 2 == 0) {
                left++;
            }
            while (left < right and nums[right] % 2 == 1) {
                right--;
            }
            if (left < right) {
                swap(nums[left++], nums[right--]);
            }
        }
        return nums;
    }
}

void sortArrayByParity_test() {
    vector<int> nums, ans;
    nums = {3, 1, 2, 4};
    ans = sortArrayByParity::sortArrayByParity(nums);
    print_vector(ans);
    nums = {0};
    ans = sortArrayByParity::sortArrayByParity(nums);
    print_vector(ans);
}

namespace sumSubarrayMins {
    int sumSubarrayMins(vector<int> &arr) {
        int n = arr.size();
        long long ans = 0;
        long long mod = 1e9 + 7;
        stack<int> monoStack;
        vector<int> dp(n);
        for (int i = 0; i < n; i++) {
            while (!monoStack.empty() && arr[monoStack.top()] > arr[i]) {
                monoStack.pop();
            }
            int k = monoStack.empty() ? (i + 1) : (i - monoStack.top());
            dp[i] = k * arr[i] + (monoStack.empty() ? 0 : dp[i - k]);
            ans = (ans + dp[i]) % mod;
            monoStack.emplace(i);
        }
        return ans;
    }
}

void sumSubarrayMins_test() {
    vector<int> arr;
    arr = {3, 1, 2, 4};
    cout << sumSubarrayMins::sumSubarrayMins(arr) << endl;
    arr = {11, 81, 94, 43, 3};
    cout << sumSubarrayMins::sumSubarrayMins(arr) << endl;
}

namespace smallestRangeI {
    int smallestRangeI(vector<int> &nums, int k) {
        int min_num = *min_element(nums.begin(), nums.end());
        int max_num = *max_element(nums.begin(), nums.end());
        return max_num - min_num <= 2 * k ? 0 : max_num - min_num - 2 * k;
    }
}

void smallestRangleI_test() {
    vector<int> nums;
    int k;
    nums = {1};
    k = 0;
    cout << smallestRangeI::smallestRangeI(nums, k) << endl;
    nums = {0, 10};
    k = 2;
    cout << smallestRangeI::smallestRangeI(nums, k) << endl;
    nums = {1, 3, 6};
    k = 3;
    cout << smallestRangeI::smallestRangeI(nums, k) << endl;
}

namespace smallestRangeII {
    int smallestRangeII(vector<int> &nums, int k) {
        sort(nums.begin(), nums.end());
        int mi = nums[0], ma = nums.back();
        int res = ma - mi;
        for (int i = 0; i < nums.size() - 1; i++) {
            int a = nums[i], b = nums[i + 1];
            res = min(res, max(ma - k, a + k) - min(mi + k, b - k));
        }
        return res;
    }
}

void smallestRangeII_test() {
    vector<int> nums;
    int k;
    nums = {1};
    k = 0;
    cout << smallestRangeII::smallestRangeII(nums, k) << endl;
    nums = {0, 10};
    k = 6;
    cout << smallestRangeII::smallestRangeII(nums, k) << endl;
    nums = {1, 3, 6};
    k = 3;
    cout << smallestRangeII::smallestRangeII(nums, k) << endl;
}

namespace sortArray {
    vector<int> tmp;

    void mergeSort(vector<int> &nums, int l, int r) {
        if (l >= r) return;
        int mid = (l + r) >> 1;
        mergeSort(nums, l, mid);
        mergeSort(nums, mid + 1, r);
        int i = l, j = mid + 1;
        int cnt = 0;
        while (i <= mid && j <= r) {
            if (nums[i] <= nums[j]) {
                tmp[cnt++] = nums[i++];
            } else {
                tmp[cnt++] = nums[j++];
            }
        }
        while (i <= mid) {
            tmp[cnt++] = nums[i++];
        }
        while (j <= r) {
            tmp[cnt++] = nums[j++];
        }
        for (int i = 0; i < r - l + 1; ++i) {
            nums[i + l] = tmp[i];
        }
    }

    vector<int> sortArray(vector<int> &nums) {
        tmp.resize((int) nums.size(), 0);
        mergeSort(nums, 0, (int) nums.size() - 1);
        return nums;
    }
}

void sortArray_test() {
    vector<int> nums, ans;
    nums = {5, 2, 3, 1};
    ans = sortArray::sortArray(nums);
    print_vector(ans);
    nums = {5, 1, 1, 2, 0, 0};
    ans = sortArray::sortArray(nums);
    print_vector(ans);
}

namespace catMouseGame {
    const int MOUSE_TURN = 0, CAT_TURN = 1;
    const int DRAW = 0, MOUSE_WIN = 1, CAT_WIN = 2;
    vector<vector<int>> graph;
    vector<vector<vector<int>>> degrees;
    vector<vector<vector<int>>> results;

    vector<tuple<int, int, int>> GetPrevStates(int mouse, int cat, int turn) {
        vector<tuple<int, int, int>> prevStates;
        int prevTurn = turn == MOUSE_TURN ? CAT_TURN : MOUSE_TURN;
        if (prevTurn == MOUSE_TURN) {
            for (int &prev : graph[mouse]) {
                prevStates.emplace_back(prev, cat, prevTurn);
            }
        } else {
            for (int &prev : graph[cat]) {
                if (prev != 0) {
                    prevStates.emplace_back(mouse, prev, prevTurn);
                }
            }
        }
        return prevStates;
    }

    int catMouseGame(vector<vector<int>> &graph) {
        int n = graph.size();
        graph = graph;
        degrees = vector<vector<vector<int>>>(n, vector<vector<int>>(n, vector<int>(2)));
        results = vector<vector<vector<int>>>(n, vector<vector<int>>(n, vector<int>(2)));
        queue<tuple<int, int, int>> qu;

        for (int i = 0; i < n; i++) {
            for (int j = 1; j < n; j++) {
                degrees[i][j][MOUSE_TURN] = graph[i].size();
                degrees[i][j][CAT_TURN] = graph[j].size();
            }
        }
        for (int node : graph[0]) {
            for (int i = 0; i < n; i++) {
                degrees[i][node][CAT_TURN]--;
            }
        }
        for (int j = 1; j < n; j++) {
            results[0][j][MOUSE_TURN] = MOUSE_WIN;
            results[0][j][CAT_TURN] = MOUSE_WIN;
            qu.emplace(0, j, MOUSE_TURN);
            qu.emplace(0, j, CAT_TURN);
        }
        for (int i = 1; i < n; i++) {
            results[i][i][MOUSE_TURN] = CAT_WIN;
            results[i][i][CAT_TURN] = CAT_WIN;
            qu.emplace(i, i, MOUSE_TURN);
            qu.emplace(i, i, CAT_TURN);
        }
        while (!qu.empty()) {
            auto[mouse, cat, turn] = qu.front();
            qu.pop();
            int result = results[mouse][cat][turn];
            vector<tuple<int, int, int>> prevStates = GetPrevStates(mouse, cat, turn);
            for (auto &[prevMouse, prevCat, prevTurn] : prevStates) {
                if (results[prevMouse][prevCat][prevTurn] == DRAW) {
                    bool canWin = (result == MOUSE_WIN && prevTurn == MOUSE_TURN) ||
                                  (result == CAT_WIN && prevTurn == CAT_TURN);
                    if (canWin) {
                        results[prevMouse][prevCat][prevTurn] = result;
                        qu.emplace(prevMouse, prevCat, prevTurn);
                    } else if (--degrees[prevMouse][prevCat][prevTurn] == 0) {
                        int loseResult = prevTurn == MOUSE_TURN ? CAT_WIN : MOUSE_WIN;
                        results[prevMouse][prevCat][prevTurn] = loseResult;
                        qu.emplace(prevMouse, prevCat, prevTurn);
                    }
                }
            }
        }
        return results[1][2][MOUSE_TURN];
    }

}

void catMouseGame_test() {
    vector<vector<int>> ga;
    ga.push_back({2, 5});
    ga.push_back({3});
    ga.push_back({0, 4, 5});
    ga.push_back({1, 4, 5});
    ga.push_back({2, 3});
    ga.push_back({0, 2, 3});

    cout << catMouseGame::catMouseGame(ga) << endl;
    ga = {{1, 3},
          {0},
          {3},
          {0, 2}};
    cout << catMouseGame::catMouseGame(ga) << endl;
}

namespace hasGroupsSizeX {
    int cnt[10000];

    bool hasGroupsSizeX(vector<int> &deck) {
        for (auto x: deck) cnt[x]++;
        int g = -1;
        for (int i = 0; i < 10000; ++i) {
            if (cnt[i]) {
                if (~g) {
                    g = gcd(g, cnt[i]);
                } else {
                    g = cnt[i];
                }
            }
        }
        return g >= 2;
    }
}

void hasGroupsSizeX_test() {
    vector<int> deck;
    deck = {1, 2, 3, 4, 4, 3, 2, 1};
    cout << hasGroupsSizeX::hasGroupsSizeX(deck) << endl;
    deck = {1, 1, 1, 2, 2, 2, 3, 3};
    cout << hasGroupsSizeX::hasGroupsSizeX(deck) << endl;
}

namespace partitionDisjoint {
    int partitionDisjoint(vector<int> &nums) {
        int n = nums.size();
        int leftMax = nums[0], leftPos = 0, curMax = nums[0];
        for (int i = 1; i < n - 1; i++) {
            curMax = max(curMax, nums[i]);
            if (nums[i] < leftMax) {
                leftMax = curMax;
                leftPos = i;
            }
        }
        return leftPos + 1;
    }
}

void partitionDisjoint_test() {
    vector<int> nums;
    nums = {5, 0, 3, 8, 6};
    cout << partitionDisjoint::partitionDisjoint(nums) << endl;
    nums = {1, 1, 1, 0, 6, 12};
    cout << partitionDisjoint::partitionDisjoint(nums) << endl;
}

namespace wordSubsets {
    void statWord(std::string &word, std::vector<int> &cnt) {
        std::vector<int> tmp_cnt(26, 0);
        for (auto &ch : word) {
            tmp_cnt[ch - 'a']++;
        }
        for (int i = 0; i < 26; i++) {
            cnt[i] = std::max(cnt[i], tmp_cnt[i]);
        }
    }

    bool isContained(vector<int> &lhs, vector<int> &rhs) {
        for (int i = 0; i < lhs.size(); i++) {
            if (rhs[i]) {
                if (lhs[i] >= rhs[i]) {
                    continue;
                }
                return false;
            }
        }
        return true;
    }

    vector<string> wordSubsets(vector<string> &words1, vector<string> &words2) {
        vector<string> ans;
        std::vector<int> ch_cnt(26, 0);
        for (auto &word : words2) {
            statWord(word, ch_cnt);
        }

        for (auto &word : words1) {
            std::vector<int> tmp_cnt(26, 0);
            statWord(word, tmp_cnt);
            if (isContained(tmp_cnt, ch_cnt)) {
                ans.push_back(word);
            }
        }
        return ans;
    }
}

void wordSubsets_test() {
    vector<string> words1, words2, ans;
    words1 = {"amazon", "apple", "facebook", "google", "leetcode"};
    words2 = {"e", "o"};
    ans = wordSubsets::wordSubsets(words1, words2);
    print_vector(ans);
    words1 = {"amazon", "apple", "facebook", "google", "leetcode"};
    words2 = {"l", "e"};
    ans = wordSubsets::wordSubsets(words1, words2);
    print_vector(ans);
    words1 = {"amazon", "apple", "facebook", "google", "leetcode"};
    words2 = {"e", "oo"};
    ans = wordSubsets::wordSubsets(words1, words2);
    print_vector(ans);
    words1 = {"amazon", "apple", "facebook", "google", "leetcode"};
    words2 = {"lo", "eo"};
    ans = wordSubsets::wordSubsets(words1, words2);
    print_vector(ans);
    words1 = {"amazon", "apple", "facebook", "google", "leetcode"};
    words2 = {"ec", "oc", "ceo"};
    ans = wordSubsets::wordSubsets(words1, words2);
    print_vector(ans);
}

namespace reverseOnlyLetters {
    string reverseOnlyLetters(string s) {
        int l, r;
        l = 0;
        r = s.size() - 1;
        while (l < r) {
            if (isalpha(s[l]) && isalpha(s[r])) {
                swap(s[l], s[r]);
                l++;
                r--;
            }
            if (!isalpha(s[l])) {
                l++;
            }
            if (!isalpha(s[r])) {
                r--;
            }
        }
        return s;
    }
}

void reverseOnlyLetters_test() {
    string s, ans;
    s = "ab-cd";
    cout << "dc-ba:" << reverseOnlyLetters::reverseOnlyLetters(s) << endl;
    s = "a-bC-dEf-ghIj";
    cout << "j-Ih-gfE-dCba:" << reverseOnlyLetters::reverseOnlyLetters(s) << endl;
    s = "Test1ng-Leet=code-Q!";
    cout << "Qedo1ct-eeLg=ntse-T!:" << reverseOnlyLetters::reverseOnlyLetters(s) << endl;
}

namespace maxSubarraySumCircular {
    int maxSubarraySumCircular(vector<int> &nums) {
        int n = nums.size();
        vector<int> leftMax(n);
        // 对坐标为 0 处的元素单独处理，避免考虑子数组为空的情况
        leftMax[0] = nums[0];
        int leftSum = nums[0];
        int pre = nums[0];
        int res = nums[0];
        for (int i = 1; i < n; ++i) {
            pre = max(pre + nums[i], nums[i]);
            res = max(res, pre);
            leftSum += nums[i];
            leftMax[i] = max(leftMax[i - 1], leftSum);
        }
        // 从右到左枚举后缀，固定后缀，选择最大前缀
        int rightSum = 0;
        for (int i = n - 1; i > 0; --i) {
            rightSum += nums[i];
            res = max(res, rightSum + leftMax[i - 1]);
        }

        return res;
    }
}

void maxSubarraySumCircular_test() {
    vector<int> nums;
    nums = {1, -2, 3, -2};
    cout << maxSubarraySumCircular::maxSubarraySumCircular(nums) << endl;
    nums = {5, -3, 5};
    cout << maxSubarraySumCircular::maxSubarraySumCircular(nums) << endl;
    nums = {3, -2, 2, -3};
    cout << maxSubarraySumCircular::maxSubarraySumCircular(nums) << endl;
}

namespace numMusicPlaylists {
    int numMusicPlaylists(int n, int goal, int k) {
        int N = n, L = goal, K = k;
        const int mod = 1e9 + 7;

        vector<vector<int>> dp(L + 1, vector<int>(N + 1));
        dp[0][0] = 1;
        for (int i = 1; i <= L; ++i)
            for (int j = 1; j <= min(i, N); ++j) {
                dp[i][j] += 1ll * dp[i - 1][j - 1] * (N - j + 1) % mod;
                dp[i][j] += 1ll * dp[i - 1][j] * max(0, j - K) % mod;
                dp[i][j] %= mod;
            }

        return dp[L][N];
    }
}

void numMusicPlaylists_test() {
    int goal, k, n;
    goal = 3, n = 3, k = 1;
    cout << numMusicPlaylists::numMusicPlaylists(n, goal, k) << endl;
    goal = 3, n = 2, k = 0;
    cout << numMusicPlaylists::numMusicPlaylists(n, goal, k) << endl;
    goal = 3, n = 2, k = 1;
    cout << numMusicPlaylists::numMusicPlaylists(n, goal, k) << endl;
}

namespace minAddToMakeValid {
    int minAddToMakeValid(string s) {
        int ans = 0;
        int leftCount = 0;
        for (auto &c : s) {
            if (c == '(') {
                leftCount++;
            } else {
                if (leftCount > 0) {
                    leftCount--;
                } else {
                    ans++;
                }
            }
        }
        ans += leftCount;
        return ans;
    }
}

void minAddToMakeValid_test() {
    string s;
    s = "())";
    cout << minAddToMakeValid::minAddToMakeValid(s) << endl;
    s = "(((";
    cout << minAddToMakeValid::minAddToMakeValid(s) << endl;
}

namespace sortArrayByParityII {
    vector<int> sortArrayByParityII(vector<int> &nums) {
        vector<int> ans(nums.size());
        int index = 0;
        for (int i = 0; i < nums.size(); ++i) {
            if (nums[i] % 2 == 0) {
                ans[index] = nums[i];
                index += 2;
            }
        }
        index = 1;
        for (int i = 0; i < nums.size(); ++i) {
            if (nums[i] % 2 == 1) {
                ans[index] = nums[i];
                index += 2;
            }
        }
        return ans;
    }
}

void sortArrayByParityII_test() {
    vector<int> nums, ans;
    nums = {4, 2, 5, 7};
    ans = sortArrayByParityII::sortArrayByParityII(nums);
    print_vector(ans);
    nums = {2, 3};
    ans = sortArrayByParityII::sortArrayByParityII(nums);
    print_vector(ans);
}

namespace threeSumMulti {
    int threeSumMulti(vector<int> &arr, int target) {
        int mod = 1e9 + 7;  // 为防止溢出
        int n = arr.size();
        sort(arr.begin(), arr.end());
        int res = 0;
        for (int i = 0; i < n - 2; i++) {
            if (arr[i] + arr[i + 1] + arr[i + 2] > target) return res;
            if (arr[i] + arr[n - 1] + arr[n - 2] < target) continue;
            int left = i + 1, right = n - 1;
            while (left < right) {
                int sum = arr[i] + arr[left] + arr[right];
                if (sum < target) left++;
                else if (sum > target) right--;
                else if (arr[left] == arr[right]) {
                    int diff = right - left + 1;
                    res = (res + diff * (diff - 1) / 2) % mod;
                    break;
                } else {
                    int l_cnt = 1, r_cnt = 1;
                    while (arr[left + 1] == arr[left++]) l_cnt++;
                    while (arr[right - 1] == arr[right--]) r_cnt++;
                    res = (res + l_cnt * r_cnt) % mod;
                }
            }
        }
        return res;
    }
}

void threeSumMulti_test() {
    vector<int> arr;
    int target = 8;
    arr = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5};
    cout << threeSumMulti::threeSumMulti(arr, target) << endl;
    target = 5;
    arr = {1, 1, 2, 2, 2, 2};
    cout << threeSumMulti::threeSumMulti(arr, target) << endl;
}

namespace minFlipsMonoIncr {
    int minFlipsMonoIncr(string s) {
        int dp0 = 0, dp1 = 0;
        for (char c:s) {
            int dp0New = dp0, dp1New = min(dp0, dp1);
            if (c == '1') {
                dp0New++;
            } else {
                dp1New++;
            }
            dp0 = dp0New;
            dp1 = dp1New;
        }
        return min(dp0, dp1);
    }
}

void minFlipsMonoIncr_test() {
    string s;
    s = "00110";
    cout << minFlipsMonoIncr::minFlipsMonoIncr(s) << endl;
    s = "010110";
    cout << minFlipsMonoIncr::minFlipsMonoIncr(s) << endl;
    s = "00011000";
    cout << minFlipsMonoIncr::minFlipsMonoIncr(s) << endl;
}

namespace threeEqualParts {
    vector<int> threeEqualParts(vector<int> &arr) {
        int sum = accumulate(arr.begin(), arr.end(), 0);
        if (sum % 3 != 0) {
            return {-1, -1};
        }
        if (sum == 0) {
            return {0, 2};
        }

        int partial = sum / 3;
        int first = 0, second = 0, third = 0, cur = 0;
        for (int i = 0; i < arr.size(); i++) {
            if (arr[i] == 1) {
                if (cur == 0) {
                    first = i;
                } else if (cur == partial) {
                    second = i;
                } else if (cur == 2 * partial) {
                    third = i;
                }
                cur++;
            }
        }

        int len = (int) arr.size() - third;
        if (first + len <= second && second + len <= third) {
            int i = 0;
            while (third + i < arr.size()) {
                if (arr[first + i] != arr[second + i] || arr[first + i] != arr[third + i]) {
                    return {-1, -1};
                }
                i++;
            }
            return {first + len - 1, second + len};
        }
        return {-1, -1};
    }
}

void threeEqualParts_test() {
    vector<int> arr, ans;
    arr = {1, 0, 1, 0, 1};
    ans = threeEqualParts::threeEqualParts(arr);
    print_vector(ans);
    arr = {1, 1, 0, 1, 1};
    ans = threeEqualParts::threeEqualParts(arr);
    print_vector(ans);
    arr = {1, 1, 0, 0, 1};
    ans = threeEqualParts::threeEqualParts(arr);
    print_vector(ans);
}

namespace minMalwareSpread {
    void dfs(vector<vector<int>> &graph, vector<int> &initialSet, vector<int> &infectedSet, int v) {
        int n = graph.size();
        for (int u = 0; u < n; u++) {
            if (graph[v][u] == 0 || initialSet[u] == 1 || infectedSet[u] == 1) {
                continue;
            }
            infectedSet[u] = 1;
            dfs(graph, initialSet, infectedSet, u);
        }
    }

    int minMalwareSpread(vector<vector<int>> &graph, vector<int> &initial) {
        int n = graph.size();
        vector<int> initialSet(n);
        for (int v : initial) {
            initialSet[v] = 1;
        }
        vector<vector<int>> infectedBy(n);
        for (int v : initial) {
            vector<int> infectedSet(n);
            dfs(graph, initialSet, infectedSet, v);
            for (int u = 0; u < n; u++) {
                if (infectedSet[u] == 1) {
                    infectedBy[u].push_back(v);
                }
            }
        }
        vector<int> count(n);
        for (int u = 0; u < n; u++) {
            if (infectedBy[u].size() == 1) {
                count[infectedBy[u][0]]++;
            }
        }
        int res = initial[0];
        for (int v : initial) {
            if (count[v] > count[res] || count[v] == count[res] && v < res) {
                res = v;
            }
        }
        return res;
    }
}

void minMalwareSpread_test() {
    vector<vector<int>> graph;
    vector<int> initial;
    graph = {{1, 1, 0},
             {1, 1, 0},
             {0, 0, 1}}, initial = {0, 1};
    cout << minMalwareSpread::minMalwareSpread(graph, initial) << endl;
    graph = {{1, 1, 0},
             {1, 1, 1},
             {0, 1, 1}}, initial = {0, 1};
    cout << minMalwareSpread::minMalwareSpread(graph, initial) << endl;
    graph = {{1, 1, 0, 0},
             {1, 1, 1, 0},
             {0, 1, 1, 1},
             {0, 0, 1, 1}}, initial = {0, 1};
    cout << minMalwareSpread::minMalwareSpread(graph, initial) << endl;
}

namespace numUniqueEmails {
    int numUniqueEmails(vector<string> &emails) {
        unordered_set<string> emailSet;
        for (auto &email: emails) {
            string local;
            for (char c: email) {
                if (c == '+' || c == '@') {
                    break;
                }
                if (c != '.') {
                    local += c;
                }
            }
            emailSet.emplace(local + email.substr(email.find('@')));
        }
        return emailSet.size();
    }
}

void numUniqueEmails_test() {
    vector<string> emails;
    emails = {"test.email+alex@leetcode.com", "test.e.mail+bob.cathy@leetcode.com",
              "testemail+david@lee.tcode.com"};
    cout << numUniqueEmails::numUniqueEmails(emails) << endl;
    emails = {"a@leetcode.com", "b@leetcode.com", "c@leetcode.com"};
    cout << numUniqueEmails::numUniqueEmails(emails) << endl;
}

namespace numSubarraysWithSum {
    int numSubarraysWithSum(vector<int> &nums, int goal) {
        int n = nums.size();
        int left1 = 0, left2 = 0, right = 0;
        int sum1 = 0, sum2 = 0;
        int ret = 0;
        while (right < n) {
            sum1 -= nums[right];
            while (left1 <= right && sum1 > goal) {
                sum1 -= nums[left1];
                left1++;
            }
            sum2 += nums[right];
            while (left2 <= right && sum2 >= goal) {
                sum2 -= nums[left2];
                left2++;
            }
            ret += left2 - left1;
            right++;
        }
        return ret;
    }
}

void numSubarraysWithSum_test() {
    vector<int> nums;
    int goal;
    nums = {1, 0, 1, 0, 1};
    goal = 2;
    cout << numSubarraysWithSum::numSubarraysWithSum(nums, goal) << endl;
    nums = {0, 0, 0, 0, 0};
    goal = 0;
    cout << numSubarraysWithSum::numSubarraysWithSum(nums, goal) << endl;
}

namespace beautifulArray {
    unordered_map<int, vector<int>> mp;

    vector<int> f(int N) {
        vector<int> ans(N, 0);
        int t = 0;
        if (mp.find(N) != mp.end()) {
            return mp[N];
        }
        if (N != 1) {
            for (auto x : f((N + 1) / 2)) {
                ans[t++] = 2 * x - 1;
            }
            for (auto x : f(N / 2)) {
                ans[t++] = 2 * x;
            }
        } else {
            ans[0] = 1;
        }
        mp[N] = ans;
        return ans;
    }

    vector<int> beautifulArray(int n) {
        mp.clear();
        return f(n);
    }
}

void beautifulArray_test() {
    int n;
    vector<int> ans;
    n = 4;
    ans = beautifulArray::beautifulArray(n);
    print_vector(ans);
    n = 5;
    ans = beautifulArray::beautifulArray(n);
    print_vector(ans);
}

namespace shortestBridge {
    int shortestBridge(vector<vector<int>> &grid) {
        int n = grid.size();
        vector<vector<int>> dirs = {{-1, 0},
                                    {1,  0},
                                    {0,  1},
                                    {0,  -1}};
        vector<pair<int, int>> island;
        queue<pair<int, int>> qu;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (grid[i][j] == 1) {
                    qu.emplace(i, j);
                    grid[i][j] = -1;
                    while (!qu.empty()) {
                        auto[x, y] = qu.front();
                        qu.pop();
                        island.emplace_back(x, y);
                        for (int k = 0; k < 4; k++) {
                            int nx = x + dirs[k][0];
                            int ny = y + dirs[k][1];
                            if (nx >= 0 && ny >= 0 && nx < n && ny < n && grid[nx][ny] == 1) {
                                qu.emplace(nx, ny);
                                grid[nx][ny] = -1;
                            }
                        }
                    }
                    for (auto &&[x, y] : island) {
                        qu.emplace(x, y);
                    }
                    int step = 0;
                    while (!qu.empty()) {
                        int sz = qu.size();
                        for (int i = 0; i < sz; ++i) {
                            auto[x, y] = qu.front();
                            qu.pop();
                            for (int k = 0; k < 4; ++k) {
                                int nx = x + dirs[k][0];
                                int ny = y + dirs[k][1];
                                if (nx >= 0 && ny >= 0 && nx < n && ny < n) {
                                    if (grid[nx][ny] == 0) {
                                        qu.emplace(nx, ny);
                                        grid[nx][ny] = -1;
                                    } else if (grid[nx][ny] == 1) {
                                        return step;
                                    }
                                }
                            }
                        }
                        step++;
                    }
                }
            }
        }
        return 0;
    }
}

void shortestBridge_test() {
    vector<vector<int>> grid;
    int ans;
    grid = {{0, 1},
            {1, 0}};
    ans = shortestBridge::shortestBridge(grid);
    cout << ans << endl;
    grid = {{0, 1, 0},
            {0, 0, 0},
            {0, 0, 1}};
    ans = shortestBridge::shortestBridge(grid);
    cout << ans << endl;
    grid = {{1, 1, 1, 1, 1},
            {1, 0, 0, 0, 1},
            {1, 0, 1, 0, 1},
            {1, 0, 0, 0, 1},
            {1, 1, 1, 1, 1}};
    ans = shortestBridge::shortestBridge(grid);
    cout << ans << endl;
}

namespace knightDialer {
    int mod = 1e9 + 7;

    int knightDialer(int n) {
        vector<vector<int>> moves = {
                {4, 6},
                {6, 8},
                {7, 9},
                {4, 8},
                {3, 9, 0},
                {},
                {1, 7, 0},
                {2, 6},
                {1, 3},
                {2, 4}
        };
        vector<vector<int>> d(2, vector<int>(10, 0));
        fill(d[1].begin(), d[1].end(), 1);
        for (int i = 2; i <= n; i++) {
            int x = i & 1;
            for (int j = 0; j < 10; j++) {
                d[x][j] = 0;
                for (int k : moves[j]) {
                    d[x][j] = (d[x][j] + d[x ^ 1][k]) % mod;
                }
            }
        }
        int res = 0;
        for (auto x : d[n % 2]) {
            res = (res + x) % mod;
        }
        return res;
    }
}

void knightDialer_test() {
    int n;
    n = 1;
    cout << knightDialer::knightDialer(n) << endl;
    n = 2;
    cout << knightDialer::knightDialer(n) << endl;
    n = 3131;
    cout << knightDialer::knightDialer(n) << endl;
}

namespace movesToStamp {
    vector<int> movesToStamp(string stamp, string target) {
        int m = stamp.size();
        int n = target.size();

        // 使用 vector 代替原来的数组
        vector<int> indegree(n - m + 1, m);
        vector<vector<int>> graph(n);

        // 队列，使用 vector
        vector<int> queue(n - m + 1);
        int l = 0, r = 0;

        // O(n * m)，
        // 判断位置为错误的点所影响的以i位置开头的点的错误点数进行建图或进队列
        for (int i = 0; i <= n - m; ++i) {
            // i开头....(m个)
            // i+0 i+1 i+m-1
            for (int j = 0; j < m; ++j) {
                if (target[i + j] == stamp[j]) {
                    if (--indegree[i] == 0) {
                        queue[r++] = i;
                    }
                } else {
                    // i + j
                    // from : 错误的位置
                    // to : i开头的下标
                    graph[i + j].push_back(i);
                }
            }
        }
        // 以i开头后取消的同一个位置的取消错误不要重复统计
        // 访问标记，使用 vector
        vector<bool> visited(n, false);
        vector<int> path;

        // 队列处理
        while (l < r) {
            int cur = queue[l++];
            path.push_back(cur);

            for (int i = 0; i < m; ++i) {
                // cur + i即，以i位置开头的点向后数出 为m个的
                // 它能修正的位置,并且去清算删除它的影响
                if (!visited[cur + i]) {
                    visited[cur + i] = true;
                    for (int next : graph[cur + i]) {
                        if (--indegree[next] == 0) {
                            queue[r++] = next;
                        }
                    }
                }
            }
        }

        // 如果路径大小没有达到应有的数量，返回空数组
        if (path.size() != n - m + 1) {
            return {};
        }

        // 逆序调整路径
        reverse(path.begin(), path.end());
        return path;
    }
}

void movesToStamp_test() {
    string stamp, target;
    vector<int> ans;
    stamp = "abc", target = "ababc";
    ans = movesToStamp::movesToStamp(stamp, target);
    print_vector(ans);
    stamp = "abca", target = "aabcaca";
    ans = movesToStamp::movesToStamp(stamp, target);
    print_vector(ans);
}

namespace shortestSuperstring {
    int getOverlap(string a, string b) {
        int n = a.size();
        int m = b.size();
        int t = min(n, m);
        for (int i = n; i > 0; --i) {
            string at = a.substr(n - i, i);
            string bt = b.substr(0, i);
            if (at == bt)
                return i;
        }
        return 0;
    }

    string shortestSuperstring(vector<string> &words) {
        string ret;
        int n = words.size();
        vector<vector<int>> overlap(n, vector<int>(n, 0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                overlap[i][j] = getOverlap(words[i], words[j]);
            }
            ret += words[i];
        }
        vector<vector<string>> dp(1 << n, vector<string>(n, ret));
        for (uint32_t state = 1; state < (1 << n); state++) {
            for (int j = 0; j < n; j++) {
                if (state >> j & 1) {
                    uint32_t preState = (state ^ (1 << j));
                    if (preState == 0) {
                        dp[state][j] = words[j];
                    } else {
                        for (int i = 0; i < n; i++) {
                            if (state >> i & 1) {
                                if (j == i)
                                    continue;
                                int t = dp[preState][i].size() + words[j].size() - overlap[i][j];
                                if (dp[state][j].size() > t)
                                    dp[state][j] = dp[preState][i] +
                                                   words[j].substr(overlap[i][j], words[j].size() - overlap[i][j]);
                            }
                        }
                    }
                }
            }
        }
        uint32_t endNum = (1 << n) - 1;
        for (int i = 0; i < n; i++) {
            string s = dp[endNum][i];
            if (ret.size() > s.size())
                ret = s;
        }
        return ret;
    }
}

void shortestSuperstring_test() {
    vector<string> words;
    words = {"alex", "loves", "leetcode"};
    cout << "alexlovesleetcode," << shortestSuperstring::shortestSuperstring(words) << endl;

    words = {"catg", "ctaagt", "gcta", "ttca", "atgcatc"};
    cout << "gctaagttcatgcatc," << shortestSuperstring::shortestSuperstring(words) << endl;
}

namespace minDeletionSize {
    int minDeletionSize(vector<string> &strs) {
        int row = strs.size();
        int col = strs[0].size();
        int ans = 0;
        for (int j = 0; j < col; ++j) {
            for (int i = 1; i < row; ++i) {
                if (strs[i - 1][j] > strs[i][j]) {
                    ans++;
                    break;
                }
            }
        }
        return ans;
    }
}

void minDeletionSize_test() {
    vector<string> strs;
    strs = {"cba", "daf", "ghi"};
    cout << minDeletionSize::minDeletionSize(strs) << endl;
    strs = {"a", "b"};
    cout << minDeletionSize::minDeletionSize(strs) << endl;
    strs = {"zyx", "wvu", "tsr"};
    cout << minDeletionSize::minDeletionSize(strs) << endl;
}

namespace minIncrementForUnique {
    int minIncrementForUnique(vector<int> &nums) {
        int cnt[800000] = {0};
        for (int x : nums) {
            cnt[x]++;
        }
        int ans = 0, taken = 0;
        for (int i = 0; i < 800000; ++i) {
            if (cnt[i] >= 2) {
                taken += cnt[i] - 1;
                ans -= i * (cnt[i] - 1);
            } else if (taken > 0 && cnt[i] == 0) {
                taken--;
                ans += i;
            }
        }
        return ans;
    }
}

void minIncrementForUnique_test() {
    vector<int> nums;
    nums = {1, 2, 2};
    cout << minIncrementForUnique::minIncrementForUnique(nums) << endl;
    nums = {3, 2, 1, 2, 1, 7};
    cout << minIncrementForUnique::minIncrementForUnique(nums) << endl;
}

namespace validateStackSequences {
    bool validateStackSequences(vector<int> &pushed, vector<int> &popped) {
        int pu = 0, po = 0;
        stack<int> stk;
        stk.push(pushed[pu]);
        pu++;
        while (po < popped.size()) {
            if (stk.empty() && pu == pushed.size()) {
                break;
            }
            if (stk.empty()) {
                if (pu != pushed.size()) {
                    stk.push(pushed[pu]);
                    pu++;
                }
            }
            if (stk.top() != popped[po]) {
                if (pu >= pushed.size())
                    return false;
                stk.push(pushed[pu]);
                pu++;
            } else {
                stk.pop();
                po++;
            }
        }
        return po == popped.size() ? true : false;
    }
}

void validateStackSequences_test() {
    vector<int> pushed, popped;
    pushed = {1, 2, 3, 4, 5};
    popped = {4, 5, 3, 2, 1};
    cout << validateStackSequences::validateStackSequences(pushed, popped) << endl;
    pushed = {1, 2, 3, 4, 5};
    popped = {4, 3, 5, 1, 2};
    cout << validateStackSequences::validateStackSequences(pushed, popped) << endl;
    pushed = {1, 0};
    popped = {1, 0};
    cout << validateStackSequences::validateStackSequences(pushed, popped) << endl;
}

namespace removeStones {
    void dfs(int x, vector<vector<int>> &edge, vector<int> &vis) {
        vis[x] = true;
        for (auto &y : edge[x]) {
            if (!vis[y]) {
                dfs(y, edge, vis);
            }
        }
    }

    int removeStones(vector<vector<int>> &stones) {
        int n = stones.size();
        vector<vector<int>> edge(n);
        unordered_map<int, vector<int>> rec;
        for (int i = 0; i < n; ++i) {
            rec[stones[i][0]].push_back(i);
            rec[stones[i][1] + 10001].push_back(i);
        }
        for (auto &[_, vec] : rec) {
            int k = vec.size();
            for (int i = 1; i < k; ++i) {
                edge[vec[i - 1]].push_back(vec[i]);
                edge[vec[i]].push_back(vec[i - 1]);
            }
        }
        vector<int> vis(n);
        int num = 0;
        for (int i = 0; i < n; ++i) {
            if (!vis[i]) {
                num++;
                dfs(i, edge, vis);
            }
        }
        return n - num;
    }
}

void removeStones_test() {
    vector<vector<int>> stones;
    stones = {{0, 0},
              {0, 1},
              {1, 0},
              {1, 2},
              {2, 1},
              {2, 2}};
    cout << removeStones::removeStones(stones) << endl;
    stones = {{0, 0},
              {0, 2},
              {1, 1},
              {2, 0},
              {2, 2}};
    cout << removeStones::removeStones(stones) << endl;
}

namespace bagOfTokensScore {
    int bagOfTokensScore(vector<int> &tokens, int power) {
        sort(tokens.begin(), tokens.end());
        int n = tokens.size(), ans = 0, l = 0, r = n - 1;
        while (l <= r) {
            if (power >= tokens[l]) {
                power -= tokens[l++];
                ans++;
            } else if (l < r && ans > 0) {
                ans--;
                power += tokens[r--];
            } else break;
        }
        return ans;
    }

}

void bagOfTokensScore_test() {
    vector<int> tokens;
    int power;
    tokens = {100};
    power = 50;
    //cout << bagOfTokensScore::bagOfTokensScore(tokens, power) << endl;
    tokens = {200, 100};
    power = 150;
    cout << bagOfTokensScore::bagOfTokensScore(tokens, power) << endl;
    tokens = {100, 200, 300, 400};
    power = 200;
    cout << bagOfTokensScore::bagOfTokensScore(tokens, power) << endl;
}


namespace largestTimeFromDigits {
    string largestTimeFromDigits(vector<int> &arr) {
        sort(arr.begin(), arr.end(), greater<int>());
        int h, m;
        do {
            h = arr[0] * 10 + arr[1];
            m = arr[2] * 10 + arr[3];
            if (h < 24 && m < 60) break;
        } while (next_permutation(arr.begin(), arr.end(), greater<int>()));
        if (h < 24 && m < 60) {
            char ans[6];
            sprintf(ans, "%02d:%02d", h, m);
            return string(ans);
        } else return "";
    }
}

void largestTimeFromDigits_test() {
    vector<int> arr;
    arr = {1, 2, 3, 4};
    cout << largestTimeFromDigits::largestTimeFromDigits(arr) << endl;
    arr = {5, 5, 5, 5};
    cout << largestTimeFromDigits::largestTimeFromDigits(arr) << endl;
    arr = {0, 0, 0, 0};
    cout << largestTimeFromDigits::largestTimeFromDigits(arr) << endl;
    arr = {0, 0, 1, 0};
    cout << largestTimeFromDigits::largestTimeFromDigits(arr) << endl;
}

namespace deckRevealedIncreasing {
    vector<int> deckRevealedIncreasing(vector<int> &deck) {
        sort(deck.begin(), deck.end());
        // 如何从有序数组得到原先的数组。
        // n-1逆序添加，每次添加需保证当前牌处在当前牌顶。
        // 添加前将牌底元素添加到牌顶。
        int n = deck.size();
        deque<int> dq;
        for (int i = n - 1; i >= 0; i--) {
            if (!dq.empty()) {
                dq.push_front(dq.back());
                dq.pop_back();
            }
            dq.push_front(deck[i]);
        }
        vector<int> res(dq.begin(), dq.end());
        return res;
    }
}

void deckRevealedIncreasing_test() {
    vector<int> deck, ans;
    deck = {17, 13, 11, 2, 3, 5, 7};
    ans = deckRevealedIncreasing::deckRevealedIncreasing(deck);
    print_vector(ans);
}

#include <cfloat>

namespace minAreaFreeRect {
    double minAreaFreeRect(vector<vector<int>> &points) {
        int n = points.size();
        set<pair<int, int>> points_set;
        for (auto &v : points) {
            points_set.insert(make_pair(v[0], v[1]));
        }
        double res = DBL_MAX;
        for (int i = 0; i < n; ++i) {
            int x1 = points[i][0], y1 = points[i][1];
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;
                int x2 = points[j][0], y2 = points[j][1];
                for (int k = j + 1; k < n; ++k) {
                    if (k == i)
                        continue;
                    int x3 = points[k][0], y3 = points[k][1];
                    int x4 = x2 + x3 - x1, y4 = y2 + y3 - y1;
                    if (points_set.count(pair<int, int>{x4, y4}) != 0) {  //p4存在
                        vector<int> v21{x2 - x1, y2 - y1};
                        vector<int> v31{x3 - x1, y3 - y1};
                        if (v21[0] * v31[0] + v21[1] * v31[1] == 0) {
                            double cur_area = pow(pow(v21[0], 2) + pow(v21[1], 2), 0.5) *
                                              pow(pow(v31[0], 2) + pow(v31[1], 2), 0.5);
                            if (cur_area < res)
                                res = cur_area;
                        }
                    }
                }
            }
        }
        return res != DBL_MAX ? res : 0;
    }
}

void minAreaFreeRect_test() {
    vector<vector<int>> points;
    points = {{1, 2},
              {2, 1},
              {1, 0},
              {0, 1}};
    cout << minAreaFreeRect::minAreaFreeRect(points) << endl;
    points = {{0, 1},
              {2, 1},
              {1, 1},
              {1, 0},
              {2, 0}};
    cout << minAreaFreeRect::minAreaFreeRect(points) << endl;
    points = {{0, 3},
              {1, 2},
              {3, 1},
              {1, 3},
              {2, 1}};
    cout << minAreaFreeRect::minAreaFreeRect(points) << endl;
    points = {{3, 1},
              {1, 1},
              {0, 1},
              {2, 1},
              {3, 3},
              {3, 2},
              {0, 2},
              {2, 3}};
    cout << minAreaFreeRect::minAreaFreeRect(points) << endl;
}

namespace leastOpsExpressTarget {
    // target到对应的数量映射
    unordered_map<int, int> target2num;

    // 递归计算的函数
    int dfs(int x, int target) {
        if (target2num.find(target) != target2num.end()) {
            return target2num[target];
        }
        if (x == target) {
            return 0;
        } else if (x > target) {
            return min(2 * target - 1, 2 * (x - target));
        } else {
            int p = 0;
            long xp = x;;
            while (xp < target) {
                xp *= x;
                ++p;
            }
            if (xp - target >= target) {
                return p - 1 + 1 + dfs(x, target - xp / x);
            } else {
                target2num[target] = min(p - 1 + dfs(x, target - xp / x), p + dfs(x, xp - target)) + 1;
                return target2num[target];
            }
        }

    }

    int leastOpsExpressTarget(int x, int target) {
        // 如果相等，那么就是0
        return x != target ? dfs(x, target) : 0;
    }
}

void leastOpsExpressTarget_test() {
    int x, target;
    x = 3, target = 19;
    cout << leastOpsExpressTarget::leastOpsExpressTarget(x, target) << endl;
    x = 5, target = 501;
    cout << leastOpsExpressTarget::leastOpsExpressTarget(x, target) << endl;
    x = 100, target = 100000000;
    cout << leastOpsExpressTarget::leastOpsExpressTarget(x, target) << endl;
}

namespace numsSameConsecDiff {
    string path;
    vector<int> ret;

    void dfs(int n, int index, int k) {
        // 找到了一个符合要求的答案
        if (index == n) {
            ret.push_back(stoi(path));
            return;
        }

        for (int i = 0; i <= 9; ++i) {
            // i 与前一位的差值相差 k，则 index 位置可以放入 i，然后继续往下递归
            if ((path[index - 1] - '0' - i == k) || (i + '0' - path[index - 1] == k)) {

                path += '0' + i;
                dfs(n, index + 1, k);
                path.pop_back();
            }
        }
    }

    vector<int> numsSameConsecDiff(int n, int k) {
        path = "";
        ret = {};
        for (int i = 1; i <= 9; ++i) {
            // 先确定第一位避免前导零
            path += '0' + i;
            dfs(n, 1, k);
            path.pop_back();
        }
        return ret;
    }
}

void numsSameConsecDiff_test() {
    int n, k;
    vector<int> ans;
    n = 3;
    k = 7;
    ans = numsSameConsecDiff::numsSameConsecDiff(n, k);
    print_vector(ans);
    n = 2;
    k = 1;
    ans = numsSameConsecDiff::numsSameConsecDiff(n, k);
    print_vector(ans);
    n = 2;
    k = 0;
    ans = numsSameConsecDiff::numsSameConsecDiff(n, k);
    print_vector(ans);
    n = 2;
    k = 2;
    ans = numsSameConsecDiff::numsSameConsecDiff(n, k);
    print_vector(ans);
}

namespace minCameraCover {
    struct Status {
        int a; // root 必须放置摄像头的情况下，覆盖整棵树需要的摄像头数目
        int b; // 覆盖整棵树需要的摄像头数目，无论 root 是否放置摄像头
        int c; // 覆盖两棵子树需要的摄像头数目，无论节点 root 本身是否被监控到
    };

    Status dfs(TreeNode::TreeNode *root) {
        if (!root) {
            return {INT_MAX / 2, 0, 0};
        }
        auto[la, lb, lc] = dfs(root->left);
        auto[ra, rb, rc] = dfs(root->right);
        int a = lc + rc + 1;
        int b = min(a, min(la + rb, ra + lb));
        int c = min(a, lb + rb);
        return {a, b, c};
    }

    int minCameraCover(TreeNode::TreeNode *root) {
        auto[a, b, c] = dfs(root);
        return b;
    }
}

void minCameraCover_test() {
    vector<int> vals;
    TreeNode::TreeNode *root;
    vals = {0, 0, -1, 0, 0};
    root = create_treenode(vals, true);
    cout << minCameraCover::minCameraCover(root) << endl;
    vals = {0, 0, -1, 0, -1, 0, -1, -1, 0};
    root = create_treenode(vals, true);
    cout << minCameraCover::minCameraCover(root) << endl;
}

namespace pancakeSort {
    vector<int> pancakeSort(vector<int>& arr) {
        vector<int> ret;
        for (int n = arr.size(); n > 1; n--) {
            int index = max_element(arr.begin(), arr.begin() + n) - arr.begin();
            if (index == n - 1) {
                continue;
            }
            reverse(arr.begin(), arr.begin() + index + 1);
            reverse(arr.begin(), arr.begin() + n);
            ret.push_back(index + 1);
            ret.push_back(n);
        }
        return ret;
    }
}

void pancakeSort_test() {
    vector<int>arr, ans;
    arr = {3,2,4,1};
    ans = pancakeSort::pancakeSort(arr);
    print_vector(ans);
    cout << "+++++++++++++" << endl;
    arr = {1,2,3};
    ans = pancakeSort::pancakeSort(arr);
    print_vector(ans);
    cout << "+++++++++++++" << endl;
}

namespace powerfulIntegers {
    vector<int> powerfulIntegers(int x, int y, int bound) {
        unordered_set<int> cnt;
        int value1 = 1;
        for (int i = 0; i < 21; i++) {
            int value2 = 1;
            for (int j = 0; j < 21; j++) {
                int value = value1 + value2;
                if (value <= bound) {
                    cnt.emplace(value);
                } else {
                    break;
                }
                value2 *= y;
            }
            if (value1 > bound) {
                break;
            }
            value1 *= x;
        }
        return vector<int>(cnt.begin(), cnt.end());
    }
}

void powerfulIntegers_test() {
    int x, y, bound;
    vector<int>ans;
    x = 2, y = 3, bound = 10;
    ans = powerfulIntegers::powerfulIntegers(x, y, bound);
    print_vector(ans);
    cout << "+++++++++++++" << endl;
    x = 3, y = 5, bound = 15;
    ans = powerfulIntegers::powerfulIntegers(x, y, bound);
    print_vector(ans);
    cout << "+++++++++++++" << endl;

}

namespace flipMatchVoyage {
    bool dfs(TreeNode::TreeNode* root, vector<int>&voyage, int &i, vector<int>&res) {
        if (!root) {
            return true;
        }
        if (root->val != voyage[i++]) {
            return false;
        }
        if (root->left && root->left->val != voyage[i]) {
            res.push_back(root->val);
            return dfs(root->right, voyage, i, res) && dfs(root->left, voyage, i, res);
        }
        return dfs(root->left, voyage, i, res) && dfs(root->right, voyage, i, res);
    }
    vector<int> flipMatchVoyage(TreeNode::TreeNode* root, vector<int>& voyage) {
        vector<int>res;
        int i = 0;
        if (dfs(root, voyage, i, res)) {
            return res;
        }
        return {-1};
    }
}

void flipMatchVoyage_test(){
    vector<int>vals, voyage, ans;
    TreeNode::TreeNode *root;
    vals = {1,2};
    voyage = {2,1};
    root = create_treenode(vals);
    ans = flipMatchVoyage::flipMatchVoyage(root, voyage);
    print_vector(ans);
    cout << "+++++++++" << endl;
    vals = {1,2,3};
    voyage = {1,3,2};
    root = create_treenode(vals);
    ans = flipMatchVoyage::flipMatchVoyage(root, voyage);
    print_vector(ans);
    cout << "+++++++++" << endl;
    vals = {1,2,3};
    voyage = {1,2,3};
    root = create_treenode(vals);
    ans = flipMatchVoyage::flipMatchVoyage(root, voyage);
    print_vector(ans);
    cout << "+++++++++" << endl;
}

namespace isRationalEqual{
#define x first
#define y second
    typedef unsigned long long ULL;

    typedef pair<ULL, ULL> PLL;
    ULL gcd(ULL a, ULL b){

        return b ? gcd(b, a%b) : a;

    }

    PLL simple(PLL a){//化简分数

        ULL t = gcd(a.x,a.y);

        return {a.x/t, a.y/t};

    }

    PLL add(PLL a, PLL b){//两个分数相加

        a = simple(a),b = simple(b);

        ULL down = a.y*b.y;

        ULL up = a.x*b.y + b.x*a.y;

        return simple({up,down});

    }

    PLL convert(string &s){//将小数转为分数

        PLL a = {0,1}, b = {0,1}, c = {0, 0};//整数部分 不循环小数部分 循环小数部分

        int i = 0;

        //例如25.01(52)

        while(i < s.size() && s[i]!='.') a.x = a.x*10 + s[i]-'0', i++; // 分解出整数部分 {25,1}

        i++;//跳过小数点

        while(i < s.size() && s[i]!='(') b.x = b.x*10 + s[i]-'0', b.y = b.y*10, i++; // 分解出不循环小数部分 {1,100}

        i++;//跳过左括号

        while(i < s.size() && s[i]!=')') c.x = c.x*10 + s[i]-'0', c.y = c.y*10+9, i++;// 分解出循环小数部分 {52,99}

        c.y *= b.y;//把循环小数前面的0也计算在内 {52,9900}

        if(c.y==0) c.y = 1;//注意 可能无循环部分 分母不能为零 这里把分母设为1即可

        // cout<<a.x<<" "<<a.y<<"--"<<b.x<<" "<<b.y<<"--"<<c.x<<" "<<c.y<<endl;

        return add(add(a,b),c);

    }

    bool isRationalEqual(string s, string t) {
        auto t1 = convert(s), t2 = convert(t);
        return t1.x == t2.x && t1.y == t2.y;
    }
}

void isRationalEqual_test(){
    string s, t;
    s = "0.(52)", t = "0.5(25)";
    cout << isRationalEqual::isRationalEqual(s, t) << endl;
    s = "0.1666(6)", t = "0.166(66)";
    cout << isRationalEqual::isRationalEqual(s, t) << endl;
    s = "0.9(9)", t = "1.";
    cout << isRationalEqual::isRationalEqual(s, t) << endl;
}

namespace kClosest {
    mt19937 gen{random_device{}()};
    void random_select(vector<vector<int>>& points, int left, int right, int k) {
        int pivot_id = uniform_int_distribution<int>{left, right}(gen);
        int pivot = points[pivot_id][0] * points[pivot_id][0] + points[pivot_id][1] * points[pivot_id][1];
        swap(points[right], points[pivot_id]);
        int i = left - 1;
        for (int j = left; j < right; ++j) {
            int dist = points[j][0] * points[j][0] + points[j][1] * points[j][1];
            if (dist <= pivot) {
                ++i;
                swap(points[i], points[j]);
            }
        }
        ++i;
        swap(points[i], points[right]);
        // [left, i-1] 都小于等于 pivot, [i+1, right] 都大于 pivot
        if (k < i - left + 1) {
            random_select(points, left, i - 1, k);
        }
        else if (k > i - left + 1) {
            random_select(points, i + 1, right, k - (i - left + 1));
        }
    }

    vector<vector<int>> kClosest(vector<vector<int>>& points, int k) {
        int n = points.size();
        random_select(points, 0, n - 1, k);
        return {points.begin(), points.begin() + k};
    }
}

void kClosest_test(){
    vector<vector<int>> points, ans;
    int k;
    points = {{1,3},{-2,2}};
    k = 1;
    ans = kClosest::kClosest(points, k);
    for (auto p : ans) {
        cout << p[0] << "," << p[1]<<";";
    }
    cout << endl;
    cout << "++++++++++++++++" << endl;
    points = {{3,3},{5,-1},{-2,4}};
    k = 2;
    ans = kClosest::kClosest(points, k);
    for (auto p : ans) {
        cout << p[0] << "," << p[1]<<";";
    }
    cout << endl;
    cout << "++++++++++++++++" << endl;
}

namespace subarraysDivByK {
    int subarraysDivByK(vector<int>& nums, int k) {
        unordered_map<int, int> record = {{0,1}};
        int sum = 0, ans = 0;
        for (int elem : nums) {
            sum += elem;
            int modules = (sum % k + k) % k;
            if (record.count(modules)) {
                ans += record[modules];
            }
            ++record[modules];
        }
        return ans;
    }
}

void subarraysDivByK_test() {
    vector<int>nums;
    int k;
    nums = {4,5,0,-2,-3,1}, k = 5;
    cout << subarraysDivByK::subarraysDivByK(nums, k) << endl;
    nums = {5}, k = 9;
    cout << subarraysDivByK::subarraysDivByK(nums, k) << endl;
}

namespace oddEvenJumps {
    vector<int> make(vector<int> & idxArr) {
        vector<int> rs(idxArr.size(), -1);
        stack<int> s;
        for(int idx: idxArr) {
            while (!s.empty() && s.top() < idx) {
                int tmp = s.top(); s.pop();
                rs[tmp] = idx;
            }
            s.push(idx);
        }
        return rs;
    }

    int oddEvenJumps(vector<int>& arr) {
        int n = arr.size();
        vector<int> iArr(n), odd(n), even(n), oddNext, evenNext;
        odd[n-1] = even[n-1] = true;

        for(int i = 0; i < n; i++) iArr[i] = i;
        sort(iArr.begin(), iArr.end(), [&](auto a, auto b){return arr[a] == arr[b] ? a < b: arr[a] < arr[b];});
        oddNext = make(iArr);
        sort(iArr.begin(), iArr.end(), [&](auto a, auto b){return arr[a] == arr[b] ? a < b: arr[a] > arr[b];});
        evenNext = make(iArr);

        for(int i = n-2; i >= 0; i--) {
            if (oddNext[i] != -1) {
                int nxt = oddNext[i];
                odd[i] = even[nxt];
            }
            if (evenNext[i] != -1) {
                int nxt = evenNext[i];
                even[i] = odd[nxt];
            }
        }
        int ans = 0;
        for (int i = 0; i < n; i++) {
            if (odd[i] == true) ans++;
        }
        return ans;
    }
}

void oddEvenJumps_test(){
    vector<int>arr;
    arr = {10,13,12,14,15};
    cout << oddEvenJumps::oddEvenJumps(arr) << endl;
    arr = {2,3,1,1,4};
    cout << oddEvenJumps::oddEvenJumps(arr) << endl;
}

namespace largestPerimeter {
    int largestPerimeter(vector<int>& nums) {
        sort(nums.begin(), nums.end());
        for (int i = (int)nums.size() - 1; i >= 2; --i){
            if (nums[i - 2] + nums[i - 1] > nums[i]) {
                return nums[i - 2] + nums[i - 1] + nums[i];
            }
        }
        return 0;
    }
}

void largestPerimeter_test() {
    vector<int>nums;
    nums = {2,1,2};
    cout << largestPerimeter::largestPerimeter(nums) << endl;
    nums = {1,2,1,10};
    cout << largestPerimeter::largestPerimeter(nums) << endl;
}

namespace sortedSquares {
    vector<int> sortedSquares(vector<int>& nums) {
        int n = nums.size();
        vector<int> ans(n);
        for (int i = 0, j = n - 1, pos = n - 1; i <= j;) {
            if (nums[i] * nums[i] > nums[j] * nums[j]) {
                ans[pos] = nums[i] * nums[i];
                ++i;
            } else {
                ans[pos] = nums[j] * nums[j];
                --j;
            }
            --pos;
        }
        return ans;
    }
}

void sortedSquares_test(){
    vector<int>nums, ans;
    nums = {-4,-1,0,3,10};
    ans = sortedSquares::sortedSquares(nums);
    print_vector(nums);
    nums = {-7,-3,2,3,11};
    ans = sortedSquares::sortedSquares(nums);
    print_vector(nums);
}

namespace maxTurbulenceSize {
    int maxTurbulenceSize(vector<int>& arr) {
        int n = arr.size();
        int ret = 1;
        int left = 0, right = 0;
        while (right < n - 1) {
            if (left == right) {
                if (arr[left] == arr[left + 1]) {
                    left++;
                }
                right++;
            } else {
                if (arr[right - 1] < arr[right] && arr[right] > arr[right + 1]) {
                    right++;
                } else if (arr[right - 1] > arr[right] && arr[right] < arr[right + 1]) {
                    right++;
                } else {
                    left = right;
                }
            }
            ret = max(ret, right - left + 1);
        }
        return ret;
    }
}

void maxTurbulenceSize_test(){
    vector<int>arr;
    arr = {9,4,2,10,7,8,8,1,9};
    cout << maxTurbulenceSize::maxTurbulenceSize(arr) << endl;
    arr = {4,8,12,16};
    cout << maxTurbulenceSize::maxTurbulenceSize(arr) << endl;
}

namespace distributeCoins {
    int distributeCoins(TreeNode::TreeNode* root) {
        int move = 0;

        function<int(const TreeNode::TreeNode *)> dfs = [&](const TreeNode::TreeNode *root) -> int {
            int moveleft = 0;
            int moveright = 0;
            if (root == nullptr) {
                return 0;
            }
            if (root->left) {
                moveleft = dfs(root->left);
            }
            if (root->right) {
                moveright = dfs(root->right);
            }
            move += abs(moveleft) + abs(moveright);
            return moveleft + moveright + root->val - 1;
        };

        dfs(root);
        return move;
    }
}

void distributeCoins_test(){
    vector<int>vals;
    TreeNode::TreeNode* root;
    vals = {3,0,0};
    root = create_treenode(vals, true);
    cout << distributeCoins::distributeCoins(root) << endl;
    vals = {0,3,0};
    root = create_treenode(vals, true);
    cout << distributeCoins::distributeCoins(root) << endl;
}

namespace countTriplets {
    int countTriplets(vector<int>& nums) {
        vector<int> cnt(1 << 16);
        for (int x: nums) {
            for (int y: nums) {
                ++cnt[x & y];
            }
        }
        int ans = 0;
        for (int x: nums) {
            for (int mask = 0; mask < (1 << 16); ++mask) {
                if ((x & mask) == 0) {
                    ans += cnt[mask];
                }
            }
        }
        return ans;
    }
}

void countTriplets_test(){
    vector<int>nums;
    nums = {2,1,3};
    cout << countTriplets::countTriplets(nums) << endl;
    nums = {0,0,0};
    cout << countTriplets::countTriplets(nums) << endl;
}

namespace mincostTickets {
    unordered_set<int> dayset;
    vector<int> costs_;
    int memo[366] = {0};

    int dp(int i) {
        if (i > 365) {
            return 0;
        }
        if (memo[i] != -1) {
            return memo[i];
        }
        if (dayset.count(i)) {
            memo[i] = min(min(dp(i + 1) + costs_[0], dp(i + 7) + costs_[1]), dp(i + 30) + costs_[2]);
        } else {
            memo[i] = dp(i + 1);
        }
        return memo[i];
    }

    int mincostTickets(vector<int>& days, vector<int>& costs) {
        costs_ = costs;
        dayset.clear();
        for (int d: days) {
            dayset.insert(d);
        }
        memset(memo, -1, sizeof(memo));
        return dp(1);
    }
}

void mincostTickets_test(){
    vector<int> days, costs;
    days = {1,4,6,7,8,20}, costs = {2,7,15};
    cout << mincostTickets::mincostTickets(days, costs) << endl;
    days = {1,2,3,4,5,6,7,8,9,10,30,31}, costs = {2,7,15};
    cout << mincostTickets::mincostTickets(days, costs) << endl;
}

namespace strWithout3a3b {
    string strWithout3a3b(int a, int b) {
        string s;
        int A = a, B = b;
        int cnta = 0, cntb = 0;
        while(A>0 || B>0)
            if (A>B) {
                if (cnta < 2) {
                    s = s + 'a';
                    cnta++;
                    A--;
                    cntb = 0;
                } else {
                    s = s + 'b';
                    cntb++;
                    B--;
                    cnta = 0;
                }
            }else {
                if (cntb < 2) {
                    s = s + 'b';
                    cntb++;
                    B--;
                    cnta = 0;
                } else {
                    s = s + 'a';
                    cnta++;
                    A--;
                    cntb = 0;
                }
            }
        return s;
    }
}

void strWithout3a3b_test() {
    int a, b;
    a = 1, b = 2;
    cout << strWithout3a3b::strWithout3a3b(a, b) << endl;
    a = 4, b = 1;
    cout << strWithout3a3b::strWithout3a3b(a, b) << endl;
}

namespace sumEvenAfterQueries {
    vector<int> sumEvenAfterQueries(vector<int>& nums, vector<vector<int>>& queries) {
        vector<int> res;
        int sum = 0;
        for(int i = 0; i < nums.size(); i++)
            if(nums[i]%2 == 0)
                sum += nums[i];
        for(int i = 0; i < queries.size(); i++)
        {
            int val = queries[i][0];
            int index = queries[i][1];
            if(nums[index] % 2 == 0)
            {
                if(val % 2 == 0)
                    sum += val;
                else
                    sum -= nums[index];
            }
            else
            {
                if(val % 2 != 0)
                    sum += nums[index]+val;
            }
            nums[index] += val;
            res.push_back(sum);
        }
        return res;
    }
}

void sumEvenAfterQueries_test(){
    vector<int>nums, ans;
    vector<vector<int>> queries;
    nums = {1,2,3,4};
    queries = {{1,0},{-3,1},{-4,0},{2,3}};
    ans  = sumEvenAfterQueries::sumEvenAfterQueries(nums, queries);
    print_vector(ans);
}

namespace intervalIntersection {
    vector<vector<int>> intervalIntersection(vector<vector<int>>& firstList, vector<vector<int>>& secondList) {
        int i = 0, j = 0;
        vector<vector<int>> res;
        while(i < firstList.size() && j < secondList.size()){
            int low = max(firstList[i][0], secondList[j][0]);
            int high = min(firstList[i][1], secondList[j][1]);
            if(low <= high){
                res.push_back({low, high});
            }
            if(firstList[i][1] < secondList[j][1]){
                i++;
            }
            else{
                j++;
            }
        }
        return res;
    }
}

void intervalIntersection_test(){
    vector<vector<int>>firstList, secondList, ans;
    firstList = firstList = {{0,2},{5,10},{13,23},{24,25}}, secondList = {{1,5},{8,12},{15,24},{25,26}};
    ans = intervalIntersection::intervalIntersection(firstList, secondList);
    for (auto list : ans) {
        cout << list[0] << " " << list[1] << ", ";
    }
    cout << endl;
    firstList = {{1,3},{5,9}}, secondList = {};
    ans = intervalIntersection::intervalIntersection(firstList, secondList);
    for (auto list : ans) {
        cout << list[0] << " " << list[1] << ", ";
    }
    cout << endl;
    firstList = {}, secondList = {{4,8},{10,12}};
    ans = intervalIntersection::intervalIntersection(firstList, secondList);
    for (auto list : ans) {
        cout << list[0] << " " << list[1] << ", ";
    }
    cout << endl;
    firstList = {{1,7}}, secondList = {{3,10}};
    ans = intervalIntersection::intervalIntersection(firstList, secondList);
    for (auto list : ans) {
        cout << list[0] << " " << list[1] << ", ";
    }
    cout << endl;
}

namespace verticalTraversal {
    vector<vector<int>> verticalTraversal(TreeNode::TreeNode* root) {
        vector<tuple<int, int, int>> nodes;
        function<void(TreeNode::TreeNode*, int, int)> dfs = [&](TreeNode::TreeNode* node, int row, int col) {
            if (node) {
                nodes.emplace_back(col, row, node->val);
                dfs(node->left, row + 1, col - 1);
                dfs(node->right, row + 1, col + 1);
            }
        };
        dfs(root, 0, 0);
        sort(nodes.begin(), nodes.end());
        vector<vector<int>> ans;
        int lastcol = INT_MIN;
        for (const auto& [col, row, value] : nodes) {
            if (col != lastcol) {
                lastcol = col;
                ans.emplace_back();
            }
            ans.back().push_back(value);
        }
        return ans;
    }
}

void verticalTraversal_test(){
    vector<int>vals;
    TreeNode::TreeNode* root;
    vals = {3,9,20,-1,-1,15,7};
    root = create_treenode(vals, true);
    vector<vector<int>> ans = verticalTraversal::verticalTraversal(root);
    for (auto list : ans) {
        print_vector(list);
    }
}

namespace smallestFromLeaf {
    string smallestFromLeaf(TreeNode::TreeNode* root) {
        string ans = "~";
        function<void(TreeNode::TreeNode*, string)> dfs = [&](TreeNode::TreeNode* node, string path) {
            if (node != nullptr) {
                path += (char)('a' + node->val);
                if (node->left == nullptr && node->right == nullptr) {
                    ans = min(ans, string(path.rbegin(), path.rend()));
                } else {
                    dfs(node->left, path);
                    dfs(node->right, path);
                }
            }
        };
        dfs(root, "");
        return ans;
    }
}

void smallestFromLeaf_test(){
    vector<int>vals;
    TreeNode::TreeNode* root;
    vals = {0,1,2,3,4,3,4};
    root = create_treenode(vals, true);
    cout << smallestFromLeaf::smallestFromLeaf(root) << endl;
    vals = {25,1,3,1,3,0,2};
    root = create_treenode(vals, true);
    cout << smallestFromLeaf::smallestFromLeaf(root) << endl;
}

namespace addToArrayForm {
    vector<int> addToArrayForm(vector<int>& nums, int k) {
        int n = nums.size();
        vector<int> ans;
        for (int i = n - 1; i >= 0 || k > 0; --i) {
            if (i >= 0) {
                k += nums[i];
            }
            ans.push_back(k % 10);
            k /= 10;
        }
        reverse(ans.begin(), ans.end());
        return ans;
    }
 }

void addToArrayForm_test(){
    vector<int>nums, ans;
    int k;
    nums = {9,9,9,9,9,9,9,9,9,9};
    k = 1;
    ans = addToArrayForm::addToArrayForm(nums, k);
    print_vector(ans);
    nums = {1,2,0,0};
    k = 34;
    ans = addToArrayForm::addToArrayForm(nums, k);
    print_vector(ans);
    nums = {2,7,4};
    k = 181;
    ans = addToArrayForm::addToArrayForm(nums, k);
    print_vector(ans);
    nums = {2,1,5};
    k = 806;
    ans = addToArrayForm::addToArrayForm(nums, k);
    print_vector(ans);
}

namespace equationsPossible {
    class UnionFind {
    private:
        vector<int> parent;

    public:
        UnionFind() {
            parent.resize(26);
            iota(parent.begin(), parent.end(), 0);
        }

        int find(int index) {
            if (index == parent[index]) {
                return index;
            }
            parent[index] = find(parent[index]);
            return parent[index];
        }

        void unite(int index1, int index2) {
            parent[find(index1)] = find(index2);
        }
    };

    bool equationsPossible(vector<string>& equations) {
        UnionFind uf;
        for (const string& str: equations) {
            if (str[1] == '=') {
                int index1 = str[0] - 'a';
                int index2 = str[3] - 'a';
                uf.unite(index1, index2);
            }
        }
        for (const string& str: equations) {
            if (str[1] == '!') {
                int index1 = str[0] - 'a';
                int index2 = str[3] - 'a';
                if (uf.find(index1) == uf.find(index2)) {
                    return false;
                }
            }
        }
        return true;
    }
}

void equationsPossible_test(){
    vector<string> equations;
    bool ans;
    equations = {"a==b","b!=a"};
    ans = equationsPossible::equationsPossible(equations);
    cout << ans << endl;
    equations = {"b==a","a==b"};
    ans = equationsPossible::equationsPossible(equations);
    cout << ans << endl;
    equations = {"c==c","b==d","x!=z"};
    ans = equationsPossible::equationsPossible(equations);
    cout << ans << endl;
    equations = {"c==c","b==c","c==a"};
    ans = equationsPossible::equationsPossible(equations);
    cout << ans << endl;
    equations = {"a==b", "b!=c", "c==a"};
    ans = equationsPossible::equationsPossible(equations);
    cout << ans << endl;
}

namespace brokenCalc {
    int brokenCalc(int startValue, int target) {
        if (startValue >= target) {
            return startValue - target;
        }
        int count = 0;
        while (startValue != target) {

            if (target > startValue) {
                if (target % 2 == 1) {
                    target = target + 1;
                } else {
                    target = target / 2;
                }
                count++;
            } else if (target < startValue) {
                count = count + startValue - target;
                target = startValue;
            }
        }

        return count;
    }
}

void brokenCalc_test(){
    int startValue, target;
    startValue = 2, target = 3;
    cout << brokenCalc::brokenCalc(startValue, target) << endl;
    startValue = 5, target = 8;
    cout << brokenCalc::brokenCalc(startValue, target) << endl;
    startValue = 3, target = 10;
    cout << brokenCalc::brokenCalc(startValue, target) << endl;
}

namespace subarraysWithKDistinct {
    int getMostDistinct(vector<int>& nums, int k) {
        unordered_map<int, int> mp;
        int left = 0, right = 0, ret = 0;
        while (right < nums.size()) {
            ++mp[nums[right++]];
            while (mp.size() > k) {
                --mp[nums[left]];
                if (mp[nums[left]] == 0) mp.erase(nums[left]);
                ++left;
            }
            // 如果这里改成 ret = max(ret, right - left)，那么此函数就是 LeetCode 904 题的解：求长度最大的子数组（此子数组中包含不同整数个数最多为K）
            ret += right - left;
        }
        return ret;
    }
    int subarraysWithKDistinct(vector<int>& nums, int k) {
        return getMostDistinct(nums, k) - getMostDistinct(nums, k - 1);
    }
}

void subarraysWithKDistinct_test(){
    vector<int>nums;
    int k;
    nums = {1,2,1,2,3};
    k = 2;
    cout << subarraysWithKDistinct::subarraysWithKDistinct(nums, k) << endl;
    nums = {1,2,1,3,4};
    k = 3;
    cout << subarraysWithKDistinct::subarraysWithKDistinct(nums, k) << endl;
}

int main() {
    subarraysWithKDistinct_test();
    {
        //brokenCalc_test();

        //equationsPossible_test();

        //addToArrayForm_test();

        //smallestFromLeaf_test();

        //verticalTraversal_test();
    
        //intervalIntersection_test();
    
        //sumEvenAfterQueries_test();

        //strWithout3a3b_test();

        //mincostTickets_test();

        //countTriplets_test();

        //distributeCoins_test();

        //maxTurbulenceSize_test();

        //sortedSquares_test();

        //largestPerimeter_test();

        //oddEvenJumps_test();

        //subarraysDivByK_test();

        //kClosest_test();

        //isRationalEqual_test();

        //flipMatchVoyage_test();

        //powerfulIntegers_test();

        //pancakeSort_test();

        //minCameraCover_test();

        //numsSameConsecDiff_test();

        //leastOpsExpressTarget_test();

        //minAreaFreeRect_test();

        //deckRevealedIncreasing_test();

        //largestTimeFromDigits_test();

        //bagOfTokensScore_test();

        //removeStones_test();

        //validateStackSequences_test();

        //minIncrementForUnique_test();

        //minDeletionSize_test();

        //shortestSuperstring_test();

        //movesToStamp_test();

        //knightDialer_test();

        //shortestBridge_test();

        //beautifulArray_test();

        //numSubarraysWithSum_test();

        //numUniqueEmails_test();

        //minMalwareSpread_test();

        //threeEqualParts_test();

        //minFlipsMonoIncr_test();

        //threeSumMulti_test();

        //sortArrayByParityII_test();

        //minAddToMakeValid_test();

        //numMusicPlaylists_test();

        //maxSubarraySumCircular_test();

        //reverseOnlyLetters_test();

        //wordSubsets_test();

        //partitionDisjoint_test();

        //hasGroupsSizeX_test();

        //catMouseGame_test();

        //sortArray_test();

        //smallestRangeII_test();

        //smallestRangleI_test();

        //sumSubarrayMins_test();

        //sortArrayByParity_test();

        //numPermsDISequence_test();

        //atMostNGivenDigitSet_test();

        //orderlyQueue_test();

        //subarrayBitwiseORs_test();

        //increasingBST_test();

        //isMonotonic_test();

        //allPossibleFBT_test();

        //numSpecialEquivGroups_test();

        //surfaceArea_test();

        //sumSubseqWidths_test();

        //findAndReplacePattern_test();

        //constructFromPrePost_test();

        //fairCandySwap_test();

        //superEggDrop_test();

        //possibleBipartition_test();

        //spiralMatrixIII_test();

        //uncommonFromSentences_test();

        //reachableNodes_test();

        //numRescueBoats_test();

        //decodeAtIndex_test();

        //peakIndexInMountainArray_test();

        //loudAndRich_test();

        //rectangleArea_test();

        //maxDistToClosest_test();

        //shortestPathLength_test();

        //isNStraightHand_test();

        //longestMountain_test();

        //backspaceCompare_test();

        //canVisitAllRooms_test();

        //numMagicSquaresInside_test();

        //numSimilarGroups_test();

        //sumOfDistancesInTree_test();

        //findReplaceString_test();

        //flipAndInvertImage_test();

        //maskPII_test();

        //largeGroupPositions_test();

        //consecutiveNumbersSum_test();

        //uniqueLetterString_test();

        //largestIsland_test();

        //maxProfitAssignment_test();

        //numFriendRequests_test();

        //flipgame_test();

        //shortestToChar_test();

        //minimumLengthEncoding_test();

        //mostCommonWord_test();

        //racecar_test();

        //numComponents_test();

        //ambiguousCoordinates_test();

        //numBusesToDestination_test();

        //pruneTree_test();

        //largestSumOfAverages_test();

        //largestTriangleArea_test();

        //xorGame_test();

        //maxIncreaseKeepingSkyline_test();

        //numberOfLines_test();

        //splitArraySameAverage_test();

        //uniqueMorseRepresentations_test();

        //minSwap_test();

        //champagneTower_test();

        //bestRotation_test();

        //allPathsSourceTarget_test();

        //rotateString_test();

        //numSubarrayBoundedMax_test();

        //validTicTacToe_test();

        //preimageSizeFZF_test();

        //numMatchingSubseq_test();

        //customSortString_test();

        //numTilings_test();

        //escapeGhosts_test();


        //rotatedDigits_test();

        //findCheapestPrice_test();

        //kthSmallestPrimeFraction_test();

        //isBipartite_test();

        //letterCasePermutation_test();

        //minDiffInBST_test();

        //numRabbits_test();

        //reachingPoints_test();

        //kthGrammar_test();

        //canTransform_test();

        //isIdealPermutation_test();

        //slidingPuzzle_test();

        //numJewelsInStones_test();

        //maxChunksToSorted_test();

        //reorganizeString_test();

        //isToeplitzMatrix_test();

        //minSwapsCouples_test();

        //orderOfLargestPlusSign_test();

        //partitionLabels_test();

        //countPrimeSetBits_test();

        //makeLargestSpecial_test();

        //intersectionSizeTwo_test();

        //reachNumber_test();

        //crackSafe_test();

        //openLock_test();

        //shortestCompletingWord_test();

        //dominantIndex_test();

        //minCostClimbingStairs_test();

        //nextGreatestLetter_test();

        //networkDelayTime_test();

        //cherryPickup_test();

        //deleteAndEarn_test();

        //dailyTemperatures_test();

        //asteroidCollision_test();

        //floodFill_test();

        //countPalindromicSubsequences_test();

        //MyCalendar_test();

        //selfDividingNumbers_test();

        //countOfAtoms_test();

        //splitListToParts_test();

        //pivotIndex_test();

        //accountsMerge_test();

        //longestWord_test();

        //smallestDistancePair_test();

        //findLength_test();

        //maxProfit_test();

        //numSubarrayProductLessThanK_test();

        //minimumDeleteSum_test();

        //binarySearch_test();

        //insertIntoBST_test();

        //searchBST_test();

        //fallingSquares_test();

        //canPartitionKSubsets_test();

        //findShortestSubArray_test();

        //countBinarySubstrings_test();

        //maxAreaOfIsland_test();

        //hasAlternatingBits_test();

        //topKFrequent692_test();

        //minStickers_test();

        //maxSumOfThreeSubarrays_test();

        //findRedundantConnection_test();

        //repeatedStringMatch_test();

        //findRedundantDirectedConnection_test();

        //calPoints_test();

        //validPalindrome_test();

        //judgePoint24_test();

        //cutOffTree_test();

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
