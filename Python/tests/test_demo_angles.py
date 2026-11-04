import ast
import unittest
from pathlib import Path

DEMOS_DIR = Path(__file__).resolve().parents[1] / "demos"


def _is_number(node, value):
    return isinstance(node, ast.Constant) and node.value == value


def _is_two_pi(node):
    return (
        isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Mult)
        and _is_number(node.left, 2)
        and isinstance(node.right, ast.Attribute)
        and isinstance(node.right.value, ast.Name)
        and node.right.value.id == "np"
        and node.right.attr == "pi"
    )


class TestDemoAngles(unittest.TestCase):
    def test_full_rotation_linspace_excludes_duplicate_endpoint(self):
        violations = []

        for demo_path in sorted(DEMOS_DIR.glob("*.py")):
            source = demo_path.read_text(encoding="utf-8")
            if "np.linspace" not in source:
                continue
            tree = ast.parse(source, filename=str(demo_path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call) or len(node.args) < 2:
                    continue
                if not (
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "np"
                    and node.func.attr == "linspace"
                    and _is_number(node.args[0], 0)
                    and _is_two_pi(node.args[1])
                ):
                    continue

                endpoint = next(
                    (keyword.value for keyword in node.keywords if keyword.arg == "endpoint"),
                    None,
                )
                if not (isinstance(endpoint, ast.Constant) and endpoint.value is False):
                    violations.append(f"{demo_path.name}:{node.lineno}")

        self.assertEqual(
            violations,
            [],
            "Full-rotation demo angles must exclude the duplicate 2*pi endpoint: "
            + ", ".join(violations),
        )


if __name__ == "__main__":
    unittest.main()
