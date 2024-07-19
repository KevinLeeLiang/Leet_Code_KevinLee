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
	if (!root) return 0;
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
	std::set_intersection(arr1.begin(), arr1.end(), arr2.begin(), arr2.end(), std::back_inserter(intersection));
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
	vector<vector<int>> param_2 = obj.getIntervals();// [[1,1]]
	for (const auto &interval : param_2) {
		cout << "[" << interval[0] << ", " << interval[1] << "] ";
	}
	cout << endl;

	obj.addNum(3);
	param_2 = obj.getIntervals();// [[1,1],[3,3]]
	for (const auto &interval : param_2) {
		cout << "[" << interval[0] << ", " << interval[1] << "] ";
	}
	cout << endl;

	obj.addNum(7);
	param_2 = obj.getIntervals();// [[1,1],[3,3],[7,7]]
	for (const auto &interval : param_2) {
		cout << "[" << interval[0] << ", " << interval[1] << "] ";
	}
	cout << endl;

	obj.addNum(2);
	param_2 = obj.getIntervals();// [[1,3],[7,7]]
	for (const auto &interval : param_2) {
		cout << "[" << interval[0] << ", " << interval[1] << "] ";
	}
	cout << endl;

	obj.addNum(6);
	param_2 = obj.getIntervals();// [[1,3],[6,7]]
	for (const auto &interval : param_2) {
		cout << "[" << interval[0] << ", " << interval[1] << "] ";
	}
	cout << endl;

	obj.addNum(8);
	param_2 = obj.getIntervals();// [[1,3],[6,8]]
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
	cout << "Maximum sum of submatrix not exceeding " << k << ": " << test.maxSumSubmatrix(matrix, k) << endl;
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
};// namespace getSum
void getSum_test() {
	int a = 5, b = 1;
	cout << "Sum of " << a << " and " << b << " is: " << getSum::getSum(a, b) << endl;
}

namespace kSmallestPairs {
	vector<vector<int>> kSmallestPairs(vector<int> &nums1, vector<int> &nums2, int k) {
		auto cmp = [&nums1, &nums2](const pair<int, int> &a, const pair<int, int> &b) {
			return nums1[a.first] + nums2[a.second] > nums1[b.first] + nums2[b.second];
		};
		priority_queue<pair<int, int>, vector<pair<int, int>>, decltype(cmp)> pq(cmp);
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
}// namespace kSmallestPairs

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
}// namespace wiggleMaxLength

void wiggleMaxLength_test() {
	//vector<int>nums = {1,17,5,10,13,15,10,5,16,8};
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
}// namespace canConstruct

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
};// namespace lexicalOrder

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

		return -1;// 如果没有唯一字符，返回 -1
	}
}// namespace firstUniqChar

void firstUniqChar_test() {
	string s = "leetcode";
	int index = firstUniqChar::firstUniqChar(s);

	if (index != -1) {
		cout << "The first unique character is '" << s[index] << "' at index " << index << "." << endl;
	} else {
		cout << "There is no unique character in the string." << endl;
	}
}

namespace lastRemaining {
	int lastRemaining(int n) {
		int a1 = 1;
		int k = 0, cnt = n, step = 1;
		while (cnt > 1) {
			if (k % 2 == 0) {// 正向
				a1 = a1 + step;
			} else {// 反向
				a1 = (cnt % 2 == 0) ? a1 : a1 + step;
			}
			k++;
			cnt = cnt >> 1;
			step = step << 1;
		}
		return a1;
	}
}// namespace lastRemaining
void lastRemaining_test() {
	int n = 9;
	cout << lastRemaining::lastRemaining(n) << endl;
}

namespace isRectangleCover {
	typedef pair<int, int> Point;

	bool isRectangleCover(vector<vector<int>> &rectangles) {
		long area = 0;
		int minX = rectangles[0][0], minY = rectangles[0][1], maxX = rectangles[0][2], maxY = rectangles[0][3];
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
		if (area != (long long) (maxX - minX) * (maxY - minY) || !cnt.count(pointMinMin) || !cnt.count(pointMinMax) ||
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
}// namespace isRectangleCover

void isRectangleCover_test() {
	vector<vector<int>> rectangles = {{1, 1, 3, 3},
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
				while (repTime--) t += o;
				stk.push_back(t);
			}
		}
		return getString(stk);
	}
}// namespace decodeString
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
}// namespace longestSubstring

void longestSubstring_test() {
	string s = "aaabb";
	int k = 3;
	cout << "字符串 " << s << " 子串中的每一字符出现次数都不少于 " << k << " 这一子串的长度 " << longestSubstring::longestSubstring(s, k) << endl;
	s = "ababbc";
	k = 2;
	cout << "字符串 " << s << " 子串中的每一字符出现次数都不少于 " << k << " 这一子串的长度 " << longestSubstring::longestSubstring(s, k) << endl;
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
}// namespace maxRotateFunction

void maxRotateFunction_test() {
	cout << atan(2.908 / 5.92) * 180 / M_PI << endl;
	vector<int> nums = {4, 3, 2, 6};
	cout << "nums 输入最大值：" << maxRotateFunction::maxRotateFunction(nums) << endl;
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
}// namespace findNthDigit

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
}// namespace removeKdigits

void removeKdigits_test() {
	string num = "1432219";
	int k = 3;
	std::cout << "num:" << num << "remove " << k << " 位数字后的最小数字是 " << removeKdigits::removeKdigits(num, k) << endl;
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
};// namespace canCross


void canCross_test() {
	std::vector<int> stones = {0, 1, 3, 5, 6, 8, 12, 17};// 石头的位置
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
			res.push_back((i % 3 == 0 ? (i % 5 == 0 ? "FizzBuzz" : "Fizz") : (i % 5 == 0 ? "Buzz" : to_string(i))));
		}
		return res;
	}
}// namespace fizzBuzz

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
		if (n < 3) return 0;
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
}// namespace numberOfArithmeticSlices

void numberOfArithmeticSlices_test() {
	vector<int> nums = {1, 2, 3, 8, 9, 10};
	nums = {1, 2, 3, 4};
	std::cout << numberOfArithmeticSlices::numberOfArithmeticSlices(nums) << std::endl;
}

namespace thirdMax {
	int thirdMax(vector<int> &nums) {
		sort(nums.begin(), nums.end(), greater<>());
		for (int i = 1, diff = 1; i < nums.size(); ++i) {
			if (nums[i] != nums[i - 1] && ++diff == 3) {// 此时 nums[i] 就是第三大的数
				return nums[i];
			}
		}
		return nums[0];
	}
}// namespace thirdMax

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
}// namespace addString
void addStrings_test() {
	string num1 = "456";
	string num2 = "77";
	std::cout << num1 << " add " << num2 << " is " << addString::addStrings(num1, num2) << std::endl;
}

namespace canPartition {
	bool canPartition(vector<int> &nums) {
		int sum = std::accumulate(nums.begin(), nums.end(), 0);
		if (sum % 2 != 0) {
			return false;// 如果数组元素和为奇数，则无法分割成两个和相等的子集
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
}// namespace canPartition

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
					if ((i == 0 || board[i - 1][j] != 'X') && (j == 0 || board[i][j - 1] != 'X')) {
						ans++;
					}
				}
			}
		}
		return ans;
	}
}// namespace countBattleships

void countBattleships_test() {
	std::vector<std::vector<char>> board = {
			{'X', '.', '.', 'X'},
			{'.', '.', '.', 'X'},
			{'.', '.', '.', 'X'}};

	std::cout << "Number of battleships: " << countBattleships::countBattleships(board) << std::endl;
}

namespace minMutation {
	int minMutation(string startGene, string endGene, vector<string> &bank) {
		unordered_set<string> bank_set(bank.begin(), bank.end());
		if (!bank_set.count(endGene)) return -1;
		queue<pair<string, int>> q;
		q.push({startGene, 0});
		char genes[] = {'A', 'C', 'G', 'T'};
		while (!q.empty()) {
			auto[current, steps] = q.front();
			q.pop();
			if (current == endGene) return steps;
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
		return -1;// 无法达到目标基因序列
	}
}// namespace minMutation

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
}// namespace countSegments

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

		Node(int _val) {
			val = _val;
		}

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
}// namespace levelOrder
void levelOrder_test() {
	levelOrder::Node *root = new levelOrder::Node(1);
	root->children.push_back(new levelOrder::Node(2));
	root->children.push_back(new levelOrder::Node(3));
	root->children.push_back(new levelOrder::Node(4));
	root->children.push_back(new levelOrder::Node(5));
	root->children[1]->children.push_back(new levelOrder::Node(6));
	root->children[1]->children.push_back(new levelOrder::Node(7));
	root->children[1]->children[1]->children.push_back(new levelOrder::Node(11));
	root->children[1]->children[1]->children[0]->children.push_back(new levelOrder::Node(14));
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

		sort(intervals.begin(), intervals.end(), [](const auto &u, const auto &v) {
			return u[1] < v[1];
		});
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
}// namespace eraseOverlapIntervals

void eraseOverlapIntervals_test() {
	vector<vector<int>> intervals = {{1, 2},
									 {2, 3},
									 {3, 4},
									 {1, 3}};
	cout << "移除 " << eraseOverlapIntervals::eraseOverlapIntervals(intervals) << " 来使剩下的区间没有重叠" << endl;
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
			auto it = lower_bound(startIntervals.begin(), startIntervals.end(), make_pair(intervals[i][1], 0));
			if (it != startIntervals.end()) {
				ans[i] = it->second;
			}
		}
		return ans;
	}
}// namespace findRightInterval

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
}// namespace pathSum

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
		k--;// 因为我们是从1开始的，所以先减去1

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
}// namespace findKthNumber
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
}// namespace arrangeCoins

void arrangeCoins_test() {
	int n = 5;
	cout << "给你一个数字 " << n << " ，计算并返回可形成完整阶梯行的总行数为：" << arrangeCoins::arrangeCoins(n) << endl;
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
}// namespace findDuplicates
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
}// namespace compress

void compress_test() {
	vector<char> chars = {'a', 'a', 'b', 'b', 'c', 'c', 'c'};
	string str1(chars.begin(), chars.end());
	cout << "chars: " << string(chars.begin(), chars.end()) << " 压缩后的数组的新长度：" << compress::compress(chars)
		 << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
	chars = {'a'};
	cout << "chars: " << string(chars.begin(), chars.end()) << " 压缩后的数组的新长度：" << compress::compress(chars)
		 << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
	chars = {'a', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b'};
	cout << "chars: " << string(chars.begin(), chars.end()) << " 压缩后的数组的新长度：" << compress::compress(chars)
		 << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
	chars = {'a', 'a', 'a', 'b', 'b', 'a', 'a'};
	cout << "chars: " << string(chars.begin(), chars.end()) << " 压缩后的数组的新长度：" << compress::compress(chars)
		 << " 压缩后的字符数组 " << string(chars.begin(), chars.end()) << endl;
	chars = {'a', 'b', 'c'};
	cout << "chars: " << string(chars.begin(), chars.end()) << " 压缩后的数组的新长度：" << compress::compress(chars)
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
};// namespace numberOfBoomerangs

void numberOfBoomerangs_test() {
	vector<vector<int>> points = {{0, 0},
								  {1, 0},
								  {2, 0}};
	cout << "points number of boomerangs is " << numberOfBoomerangs::numberOfBoomerangs(points) << endl;
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
}// namespace findDisappearedNumbers
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
}// namespace SerializingAndDeserializingForBinaryTrees

void SerializingAndDeserializingForBinaryTrees_test() {
	string token = "213";
	TreeNode::TreeNode *root = SerializingAndDeserializingForBinaryTrees::deserialize(token);
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
}// namespace deleteNode

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
}// namespace frequencySort
void frequencySort_test() {
	string s = "2a554442f544asfasssffffasss";
	cout << s << " 根据字符出现频率排序后 " << frequencySort::frequencySort(s) << endl;
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
}// namespace minMoves

void minMoves_test() {
	vector<int> nums = {1, 2, 3};
	cout << "最小操作次数：" << minMoves::minMoves(nums) << endl;
}

namespace fourSumCount {
	int fourSumCount(vector<int> &nums1, vector<int> &nums2, vector<int> &nums3, vector<int> &nums4) {
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
}// namespace fourSumCount

void fourSumCount_test() {
	vector<int> nums1 = {1, 2}, nums2 = {-2, -1}, nums3 = {-1, 2}, nums4 = {0, 2};
	cout << "四数相加为0的元组数：" << fourSumCount::fourSumCount(nums1, nums2, nums3, nums4) << endl;
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
}// namespace findContentChildren

void findContentChildren_test() {
	vector<int> g = {1, 2, 3};
	vector<int> s = {1, 1};
	cout << "有 " << findContentChildren::findContentChildren(g, s) << " 个小孩儿被满足" << endl;
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
}// namespace find132pattern

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
}// namespace hammingDistance

void hammingDistance_test() {
	int x = 1;
	int y = 4;
	cout << x << "和" << y << "的汉明距离：" << hammingDistance::hammingDistance(x, y) << endl;
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
}// namespace minMoves2

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
}// namespace islandPerimeter

void islandPerimeter_test() {
	vector<vector<int>> grid = {{0, 1, 0, 0},
								{1, 1, 1, 0},
								{0, 1, 0, 0},
								{1, 1, 0, 0}};
	cout << "岛屿grid的边界长度：" << islandPerimeter::islandPerimeter(grid) << endl;
}

namespace canIWin {
	bool canIWinHelper(int maxChoosableInteger, int desiredTotal, int chosen, std::unordered_map<int, bool> &memo) {
		if (memo.count(chosen)) {
			return memo[chosen];
		}

		for (int i = 1; i <= maxChoosableInteger; ++i) {
			int mask = 1 << (i - 1);
			if ((chosen & mask) == 0) {
				if (i >= desiredTotal || !canIWinHelper(maxChoosableInteger, desiredTotal - i, chosen | mask, memo)) {
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
}// namespace canIWin

void canIWin_test() {
	int maxChoosableInteger = 10;
	int desiredTotal = 21;
	string ans = "是";
	if (canIWin::canIWin(maxChoosableInteger, desiredTotal))
		ans = "必胜";
	else
		ans = "必败";
	cout << "可选取maxChoosableInteger：" << maxChoosableInteger << "，目标数值：" << desiredTotal << ",先手方 " << ans << endl;
}

namespace findSubstringInWraproundString {
	int findSubstringInWraproundString(string s) {
		vector<int> dp(26);
		int k = 0;
		for (int i = 0; i < s.length(); ++i) {
			if (i && (s[i] - s[i - 1] + 26) % 26 == 1) {// 字符之差为 1 或 -25
				++k;
			} else {
				k = 1;
			}
			dp[s[i] - 'a'] = max(dp[s[i] - 'a'], k);
		}
		return accumulate(dp.begin(), dp.end(), 0);
	}
}// namespace findSubstringInWraproundString

void findSubstringInWraproundString_test() {
	string s = "zab";
	cout << "字符串" << s << " 有 " << findSubstringInWraproundString::findSubstringInWraproundString(s) << " 个不同子串"
		 << endl;
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
		sort(words.begin(), words.end(), [&](const string &a, const string &b) {
			return a.size() < b.size();
		});
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
}// namespace findAllConcatenatedWordsInADict

void findAllConcatenatedWordsInADict_test() {
	vector<string> words = {"cat", "cats", "catsdogcats", "dog", "dogcatsdog", "hippopotamuses", "rat", "ratcatdogcat"};
	auto ans = findAllConcatenatedWordsInADict::findAllConcatenatedWordsInADict(words);

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
		sort(matchsticks.begin(), matchsticks.end(), greater<int>());// 减少搜索量

		vector<int> edges(4);
		return dfs(0, matchsticks, edges, totalLen / 4);
	}
}// namespace makesquare

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
		vector<vector<vector<int>>> dp(length + 1, vector<vector<int>>(m + 1, vector<int>(n + 1)));
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
}// namespace findMaxForm

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
			while (j < heaters.size() - 1 && abs(houses[i] - heaters[j]) >= abs(houses[i] - heaters[j + 1])) {
				j++;
				curDistance = min(curDistance, abs(houses[i] - heaters[j]));
			}
			ans = max(ans, curDistance);
		}
		return ans;
	}
}// namespace findRadius

void findRadius_test() {
	vector<int> houses = {1, 2, 3};
	vector<int> heaters = {2};
	//cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
	houses = {1, 2, 3, 4};
	heaters = {1, 4};
	//cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
	houses = {1, 5};
	heaters = {2};
	//cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
	houses = {1, 5};
	heaters = {10};
	cout << "最小半径，" << findRadius::findRadius(houses, heaters) << endl;
}

#include<ctime>

namespace randomlyGeneratePointsWithinACircle {

	class Solution {
		mt19937 gen{random_device{}()};
		uniform_real_distribution<double> dis;
		double xc, yc, r;

	public:
		Solution(double radius, double x_center, double y_center) : dis(-radius, radius), xc(x_center), yc(y_center),
																	r(radius) {}

		vector<double> randPoint() {
			while (true) {
				double x = dis(gen), y = dis(gen);
				if (x * x + y * y <= r * r) {
					return {xc + x, yc + y};
				}
			}
		}

	};
}

void randomlyGeneratePointsWithinACircle_test() {
	vector<double> param = {1.0, 0.0, 0.0};
	randomlyGeneratePointsWithinACircle::Solution sol(param[0], param[1], param[2]);
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
}

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
}

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
}

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
}

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
					if (i > 0 && i < curr.board.size() && curr.board[i - 1] == curr.board[i] &&
						curr.board[i] != curr.hand[j]) {
						choose = true;
					}
					if (choose) {
						string new_board = clean(curr.board.substr(0, i) + curr.hand[j] + curr.board.substr(i));
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
}

void findMinStep_test() {
	string board = "WRRBBW";
	string hand = "RB";
	cout << board << "," << hand << "," << findMinStep::findMinStep(board, hand) << endl;
	board = "WWRRBBWW";
	hand = "WRBRW";
	cout << board << "," << hand << "," << findMinStep::findMinStep(board, hand) << endl;
	board = "G";
	hand = "GGGGG";
	cout << board << "," << hand << "," << findMinStep::findMinStep(board, hand) << endl;
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
}

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
}

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
}

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
}

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

	int findMaximizedCapital(int k, int w, vector<int> &profits, vector<int> &capital) {
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
}

void findMaximizedCapital_test() {
	int k, w;
	k = 2;
	w = 0;
	vector<int> profits, captial;
	profits = {1, 2, 3};
	captial = {0, 1, 1};
	cout << findMaximizedCapital::findMaximizedCapital(k, w, profits, captial) << endl;
	k = 3;
	w = 0;
	profits = {1, 2, 3};
	captial = {0, 1, 2};
	cout << findMaximizedCapital::findMaximizedCapital(k, w, profits, captial) << endl;
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
}

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
}

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
}

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
		for (auto &[s, c]: cnt) {
			if (c == maxCnt) {
				ans.emplace_back(s);
			}
		}
		return ans;
	}
}

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
}

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
}

void longestPalindromeSubseq_test() {
	string s = "bbbab";
	cout << s << " 的最大回文子串序列个数：" << longestPalindromeSubseq::longestPalindromeSubseq(s) << endl;
	s = "cbbd";
	cout << s << " 的最大回文子串序列个数：" << longestPalindromeSubseq::longestPalindromeSubseq(s) << endl;
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
}

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
}

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
}

void findLUSlength_test() {
	string a, b;
	a = "aba", b = "cdc";
	cout << "序列a " << a << " 和序列b " << b << " 两个字符串的最长特殊序列个数：" << findLUSlength::findLUSlength(a, b) << endl;
	a = "aaa", b = "bbb";
	cout << "序列a " << a << " 和序列b " << b << " 两个字符串的最长特殊序列个数：" << findLUSlength::findLUSlength(a, b) << endl;
	a = "aaa", b = "aaa";
	cout << "序列a " << a << " 和序列b " << b << " 两个字符串的最长特殊序列个数：" << findLUSlength::findLUSlength(a, b) << endl;
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
}

void findLUSlength2_test() {
	vector<string> strs = {"aba", "cdc", "eae"};
	cout << "字符串数组：";
	print_vector(strs);
	cout << " 的最长特殊序列个数为 " << findLUSlength2::findLUSlength(strs) << endl;
	strs = {"aaa", "aaa", "aa"};
	cout << "字符串数组：";
	print_vector(strs);
	cout << " 的最长特殊序列个数为 " << findLUSlength2::findLUSlength(strs) << endl;
}

namespace checkSubarraySum {
	bool checkSubarraySum(vector<int>& nums, int k) {
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
}

void checkSubarraySum_test(){
	vector<int>nums;
	int k;
	nums = {23,2,4,6,7};
	k = 6;
	cout << "数组 " ;
	print_vector(nums);
	if (checkSubarraySum::checkSubarraySum(nums, k)) {
		cout << "存在和为 " << k << " 的连续子数组" << endl;
	} else {
		cout << "不存在和为 " << k << " 的连续子数组" << endl;
	}
	k = 6;
	cout << "数组 " ;
	print_vector(nums);
	if (checkSubarraySum::checkSubarraySum(nums, k)) {
		cout << "存在和为 " << k << " 的连续子数组" << endl;
	} else {
		cout << "不存在和为 " << k << " 的连续子数组" << endl;
	}
	k = 13;
	cout << "数组 " ;
	print_vector(nums);
	if (checkSubarraySum::checkSubarraySum(nums, k)) {
		cout << "存在和为 " << k << " 的连续子数组" << endl;
	} else {
		cout << "不存在和为 " << k << " 的连续子数组" << endl;
	}
}

int main() {
	checkSubarraySum_test();
	{
		//findLUSlength2_test();

		//findLUSlength_test();

		//change_test();

		//findMinMoves_test();

		//longestPalindromeSubseq_test();

		//findBottomLeftValue_test();

		//findFrequentTreeSum_test();

		//findRelativeRanks_test();

		//convertToBase7_test();

		//nextGreaterElements_test();

		//findMaximizedCapital_test();

		//findMode_test();

		//findWords_test();

		//findDiagonalOrder_test();

		//nextGreaterElement_test();

		//findMinStep_test();

		//predictTheWinner_test();

		//findMaxConsecutiveOnes_test();

		//medianSlidingWindow_test();

		//largestPalindrome_test();

		//randomlyGeneratePointsWithinACircle_test();

		//findRadius_test();

		//findMaxForm_test();

		//makesquare_test();

		//findAllConcatenatedWordsInADict_test();
		//findSubstringInWraproundString_test();
		//canIWin_test();
		//islandPerimeter_test();
		//minMoves2_test();
		//hammingDistance_test();
		//find132pattern_test();
		//findContentChildren_test();
		//fourSumCount_test();
		//minMoves_test();
		//frequencySort_test();
		//deleteNode_test();
		//SerializingAndDeserializingForBinaryTrees_test();
		//findDisappearedNumbers_test();
		//numberOfBoomerangs_test();
		//compress_test();
		//findDuplicates_test();
		//arrangeCoins_test();
		//findKthNumber_test();
		//pathSum_test();
		//findRightInterval_test();
		//eraseOverlapIntervals_test();
		//levelOrder_test();
		//countSegments_test();
		//minMutation_test();
		//countBattleships_test();
		//canPartition_test();
		//addStrings_test();
		//thirdMax_test();
		//rob_test();
		//countBits_test();
		//reverseVowels_test();
		//topKFrequent_test();
		//intersection_test();
		//SummaryRanges_test();
		//maxSumSubMatrix_test();
		//isPerfectSquare_test();
		//largestDivisibleSubset_test();
		//getSum_test();
		//kSmallestPairs_test();
		//canConstruct_test();
		//NestedInteger_test();
		//lexicalOrder_test();
		//firstUniqChar_test();
		//lastRemaining_test();
		//isRectangleCover_test();
		//decodeString_test();
		//longestSubstring_test();
		//maxRotateFunction_test();
		//findNthDigit_test();
		//removeKdigits_test();
		//canCross_test();
		//fizzBuzz_test();
		//numberOfArithmeticSlices_test();
	}
	return 0;
}
