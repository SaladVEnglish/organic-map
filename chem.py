# pyright: reportArgumentType=false
import re
from collections import defaultdict
from enum import Enum
from fractions import Fraction


class molecule_type(Enum):
    ATOM = 0
    RADICAL = 1
    MOLECULE = 2
    POSITIVE_ION = 3
    NEGATIVE_ION = 4


# Just a test
class valence:
    def __init__(
        self, formula, val: Fraction, cnt, type: molecule_type = molecule_type.ATOM
    ):
        self.formula = formula
        self.val = val
        self.cnt = cnt
        self.type = type

    def __int__(self) -> int:
        return int(self.val)

    def __str__(self) -> str:
        if self.val.denominator == 1:
            return f"{'-' if self.val < 0 else ''}{int_to_roman(abs(int(self.val)))}"
        return str(self.val)

    def __repr__(self) -> str:
        return f"valence({self.formula}({str(self)}) * {self.cnt})"

    def __mul__(self, other: int) -> "valence":
        return valence(self.formula, self.val, self.cnt * other, self.type)


class radical:
    def __init__(self, formula, val: Fraction, atom_valence: list[valence]):
        self.formula = formula
        self.valence = val
        self.atom_valence = atom_valence
        self.type = molecule_type.RADICAL

    def __repr__(self) -> str:
        return self.formula

    def __str__(self) -> str:
        return self.formula

    def __getitem__(self, key) -> valence:
        if isinstance(key, int):
            return (
                self.atom_valence[key - 1]
                if key
                else valence(self.formula, self.valence, 1, self.type)
            )
        for v in self.atom_valence:
            if v.formula == key:
                return v
        raise KeyError(key)

    def __len__(self) -> int:
        return len(self.atom_valence) + 1


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
    "NO3": radical(
        "NO3",
        -1,
        [
            valence("O", -2, 3, molecule_type.RADICAL),
            valence("N", +5, 1, molecule_type.RADICAL),
        ],
    ),
}

H_like_atoms = {
    "H": 1,
    "Li": 1,
    "Na": 1,
    "K": 1,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "Ca": 2,
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


def generate_chemical_formula(atom_dict):
    """
    根据原子数量字典生成化学公式
    :param atom_dict: 原子数量字典（如 {"H": 2, "O": 2}）
    :return: 化学公式字符串（如 "H2O2"）
    """
    formula = ""
    for atom, count in sorted(atom_dict.items()):
        if count == 1:
            formula += atom
        else:
            formula += f"{atom}{count}"
    return formula


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
    now_valence = 0
    involed_radicals = list[valence]()
    for i in radicals.values():
        cur_cnt = count_pattern_occurrences(formula, i.formula)
        now_valence += cur_cnt * i.valence
        if cur_cnt:
            involed_radicals.append(
                valence(i.formula, i.valence, cur_cnt, molecule_type.RADICAL)
            )
            for j in i.atom_valence:
                atom_dict[j.formula] -= cur_cnt * j.cnt
                print(j.formula, j.cnt * cur_cnt)
                if atom_dict[j.formula] <= 0:
                    del atom_dict[j.formula]

    return valence(
        generate_chemical_formula(atom_dict), -now_valence, 1, molecule_type.RADICAL
    ), involed_radicals


def simple_get_valence(formula: str, now_valence: int = 0) -> dict[str, valence]:
    atom_dict = parse_chemical_formula(formula)
    atom_cnt = len(atom_dict)
    res = dict[str, valence]()
    visit = set()
    for i in preset_valence_order.keys():
        if atom_cnt == 1:
            break
        if i in atom_dict:
            res[i] = valence(i, Fraction(preset_valence_order[i]), atom_dict[i])
            now_valence += int(atom_dict[i] * res[i].val)
            visit.add(i)
            atom_cnt -= 1
    last_atom = list(set(atom_dict.keys()) - visit)[0]
    # print(last_atom)
    res[last_atom] = valence(
        last_atom,
        -Fraction(now_valence) / Fraction(atom_dict[last_atom]),
        atom_dict[last_atom],
    )
    return res


def get_valence(formula: str, keep_radical: bool = False) -> dict[str, valence]:
    left_radical, radical_list = parse_radical(formula)
    res = simple_get_valence(left_radical.formula, -left_radical.val)
    if keep_radical:
        for i in radical_list:
            res[i.formula] = i
        return res
    else:
        for i in radical_list:
            cur_radical = radicals[i.formula]
            for j in cur_radical.atom_valence:
                if j.formula not in res:
                    res[j.formula] = j * i.cnt
                else:
                    cur_atom_cnt: int = res[j.formula].cnt + j.cnt * i.cnt
                    cur_atom_valence: Fraction = (
                        res[j.formula].cnt * res[j.formula].val + j.cnt * i.cnt * j.val
                    )
                    res[j.formula] = valence(
                        j.formula, cur_atom_valence / cur_atom_cnt, cur_atom_cnt, j.type
                    )
        return res


def int_to_roman(num):
    """
    将整数转换为罗马数字
    :param num: 整数（0-3999）
    :return: 罗马数字字符串
    """
    if num == 0:
        return "0"
    # 定义值到符号的映射（按从大到小排列）
    val_to_sym = [
        (1000, "M"),
        (900, "CM"),
        (500, "D"),
        (400, "CD"),
        (100, "C"),
        (90, "XC"),
        (50, "L"),
        (40, "XL"),
        (10, "X"),
        (9, "IX"),
        (5, "V"),
        (4, "IV"),
        (1, "I"),
    ]

    roman = ""
    for value, symbol in val_to_sym:
        # 当 num 大于等于当前值时，重复减去并添加符号
        while num >= value:
            roman += symbol
            num -= value
    return roman


class compound:
    def __init__(self, formula: str = "H2O"):
        self.formula = formula
        self.atoms = parse_chemical_formula(formula)
        self.atom_cnt = len(self.atoms)
        self.valence = get_valence(formula)

    def regen_formula(self):
        return generate_chemical_formula(self.atoms)

    def unsaturation_degree(self) -> int:
        C_cnt = 0
        H_cnt = 0
        for i in self.atoms.keys():
            if i not in H_like_atoms and i != "C" and i != "O":
                return -1
            elif i == "O":
                continue
            elif i == "C":
                C_cnt = self.atoms[i]
            else:
                H_cnt += self.atoms[i] * H_like_atoms[i]

        return (C_cnt * 2 + 2 - H_cnt) // 2

    def __repr__(self) -> str:
        return self.formula


if __name__ == "__main__":
    # print(get_valence(input()))
    # print(simple_get_valence(parse_chemical_formula(input())))
    cp = compound("H2C=CHCH3")
    print(cp.unsaturation_degree())
    print(cp.valence)
