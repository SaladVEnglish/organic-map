import re
from collections import defaultdict
from fractions import Fraction
from enum import Enum


class molecule_type(Enum):
    ATOM = 0
    RADICAL = 1
    MOLECULE = 2
    POSITIVE_ION = 3
    NEGATIVE_ION = 4


class valence:
    def __init__(self, formula, val, cnt, type: molecule_type = molecule_type.ATOM):
        self.formula = formula
        self.val = val
        self.cnt = cnt
        self.type = type


preset_valence_order = {
    "F": -1,
    "Li": 1,
    "Na": 1,
    "K": 1,
    "B": 3,
    "Al": 3,
    "H": 1,
    "O": -2,
    "Cl": -1,
    "N": -3,
    "S": -2,
    "C": 4,
    "Cu": 2,
}

changeable_valence_order = {
    "Fe": [+2, +3],
    "N": [-3, +5, +1, +2, +4],
    "S": [+6, -2, +4],
}

radicals = {
    "NO3": {
        valence("NO3", -1, 1, molecule_type.RADICAL),
        valence("O", -2, 3, molecule_type.RADICAL),
        valence("N", +5, 1, molecule_type.RADICAL),
    },
}


def parse_chemical_formula(formula):
    """
    解析化学式，返回原子及其数量的字典
    :param formula: 化学公式字符串（如 "H2O2"）
    :return: 原子数量字典（如 {"H": 2, "O": 2}）
    """
    # 使用正则表达式匹配原子、数字、括号
    # 匹配模式：原子（大写字母开头，小写字母结尾）+ 数字（可选）
    # 或者括号 + 数字（可选）
    pattern = r"([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)"
    tokens = re.findall(pattern, formula)

    stack = [defaultdict(int)]
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token[0]:  # 原子
            atom = token[0]
            count = int(token[1]) if token[1] else 1
            stack[-1][atom] += count
        elif token[2] == "(":  # 左括号
            stack.append(defaultdict(int))
        elif token[3] == ")":  # 右括号
            count = int(token[4]) if token[4] else 1
            inner = stack.pop()
            for atom, c in inner.items():
                stack[-1][atom] += c * count
        i += 1

    result = dict(stack[0])
    return result


def count_pattern_occurrences(target_str, pattern):
    """
    全局匹配目标字符串中的 pattern 变体，并累加所有匹配的数量
    变体规则：
    - pattern        → 数量+1
    - (pattern)      → 数量+1
    - (pattern)数字  → 数量+数字
    :param target_str: 待匹配的原始字符串（如 "Fe2(NO3)(NO3)4Cl"）
    :param pattern: 核心匹配模式（如 "NO3"）
    :return: 累加后的总数量（int）
    """
    # 转义pattern中的特殊字符（如()、+、*等），避免正则语法冲突
    escaped_pattern = re.escape(pattern)

    # 构建正则表达式：匹配 "pattern" / "(pattern)" / "(pattern)数字" 三种变体
    # 分组说明：
    # group1: 匹配纯pattern（无括号）
    # group2: 匹配带括号的pattern
    # group3: 匹配(pattern)后面的数字（可选）
    regex = re.compile(rf"(?:{escaped_pattern}|\({escaped_pattern}\))(\d*)")

    # 全局查找所有匹配项
    matches = regex.findall(target_str)

    total = 0
    for num_str in matches:
        # 数字为空 → 数量+1；数字不为空 → 数量+数字
        if num_str:
            total += int(num_str)
        else:
            total += 1

    return total


def parse_radical(formula):
    atom_dict = parse_chemical_formula(formula)
    for i in radicals.keys():
        pattern = r"(\()?%s((\))(\d*))?" % i


def simple_get_valence(formula: str, now_valence: int = 0) -> dict:
    atom_dict = parse_chemical_formula(formula)
    atom_cnt = len(atom_dict)
    res = dict()
    visit = set()
    for i in preset_valence_order.keys():
        if atom_cnt == 1:
            break
        if i in atom_dict:
            res[i] = Fraction(preset_valence_order[i])
            now_valence += atom_dict[i] * res[i]
            visit.add(i)
            atom_cnt -= 1
    last_atom = list(set(atom_dict.keys()) - visit)[0]
    # print(last_atom)
    res[last_atom] = -Fraction(now_valence) / Fraction(atom_dict[last_atom])
    return res


class compound:
    def __init__(self, formula: str = "H2O"):
        self.formula = formula
        self.atoms = parse_chemical_formula(formula)
        self.atom_cnt = len(self.atoms)

    def __repr__(self) -> str:
        return self.formula


if __name__ == "__main__":
    print(simple_get_valence(parse_chemical_formula(input())))
