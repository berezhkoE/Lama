#include <vector>
#include <string>
#include <variant>
#include <stack>
#include <iostream>
#include <unordered_map>
#include "../byterun/bytelib.h"
#include "../runtime/runtime.h"


namespace tkn {

    enum class Operator {
        PLUS,
        MINUS,
        MULT,
        DIV,
        MOD,
        LT,
        LE,
        GT,
        GE,
        EQ,
        NEQ,
        CONJ,
        DISJ
    };

    struct BinOp {
        Operator op;
    };

    struct Const {
        int x;
    };

    struct String {
        char *s;
    };

    struct Sexp {
        char *s;
        int n;
    };

    namespace dsg {
        struct Global {
            int x;
        };

        struct Local {
            int x;
        };

        struct Arg {
            int x;
        };

        struct Const {
            int x;
        };

        using Designation =
            std::variant<std::monostate, Global, Local, Arg, Const>;
    }

    struct Ld {
        dsg::Designation d;
    };

    struct Lda {
        dsg::Designation d;
    };

    struct St {
        dsg::Designation d;
    };

    struct Sti {
    };

    struct Sta {
    };

    struct Elem {
    };

    struct Jmp {
        int i;
    };

    struct CJmp {
        bool notZero;
        int i;
    };

    namespace call {
        struct CallRead {
        };

        struct CallWrite {
        };

        struct CallLength {
        };

        struct CallString {
        };

        struct CallArray {
            int n;
        };

        using Call = std::variant<CallRead, CallWrite, CallLength, CallString, CallArray>;
    }

    struct Begin {
        int n;
        int m;
    };

    struct End {
    };

    struct Closure {
        int n;
    };

    struct CallC {
        int x;
    };

    struct Call {
        int i;
        int n;
    };

    struct Ret {
    };

    struct Drop {
    };

    struct Dup {
    };

    struct Swap {
    };

    struct Tag {
        char *s;
        int n;
    };

    struct Array {
        int n;
    };

    struct Patt {
        std::string s;
    };

    struct Fail {
        int n;
        int m;
    };

    struct Line {
        int n;
    };

    struct CBegin {
        int n;
        int m;
    };

    using Token = std::variant<BinOp, Const, String, Sexp, Sti, Sta, Jmp, End, Ret, Drop, Dup, Swap, Elem, Ld, Lda, St, CJmp, call::Call, Begin, CBegin, Closure, CallC, Call, Tag, Array, Fail, Line, Patt>;
}

struct Scope {
    std::vector<int> args;
    std::vector<int> locals;
};

using namespace tkn;

std::vector<Token> disassemble(bytefile *bf) {
# define INT    (ip += sizeof (int), *(int*)(ip - sizeof (int)))
# define BYTE   *ip++
# define STRING get_string (bf, INT)
# define FAIL   failure ("ERROR: invalid opcode %d-%d\n", h, l)

  char *ip = bf->code_ptr;
  Operator ops[] = {Operator::PLUS, Operator::MINUS, Operator::MULT, Operator::DIV, Operator::MOD, Operator::LT,
                    Operator::LE, Operator::GT, Operator::GE, Operator::EQ, Operator::NEQ, Operator::CONJ, Operator::DISJ};
  std::string pats[] = {"=str", "#string", "#array", "#sexp", "#ref", "#val", "#fun"};
  std::string lds[] = {"LD", "LDA", "ST"};

  std::vector<Token> result;
  std::unordered_map<int, int> jmpMap;
  auto insertToJmpMap = [&]() -> void {
      jmpMap.insert({ip - bf->code_ptr - 1, result.size()});
  };

//  auto f = stdout;
  do {
    char x = BYTE,
        h = (x & 0xF0) >> 4,
        l = x & 0x0F;

    switch (h) {
      case 15:
        goto stop;

      case 0:
        //fprintf(f, "BINOP\t%d\n", ops[l - 1]);

        insertToJmpMap();
        result.emplace_back(BinOp{ops[l - 1]});
        break;

      case 1:
        switch (l) {
          case 0:
            //fprintf(f, "CONST\n");

            insertToJmpMap();
            result.emplace_back(Const{BOX(INT)});
            break;

          case 1:
            //fprintf(f, "STRING\n");

            insertToJmpMap();
            result.emplace_back(String{STRING});
            break;

          case 2:
            //fprintf(f, "SEXP\n");

            insertToJmpMap();
            result.emplace_back(Sexp{STRING, INT});
            break;

          case 3:
            //fprintf(f, "STI\n");

            insertToJmpMap();
            result.emplace_back(Sti{});
            break;

          case 4:
            //fprintf(f, "STA\n");

            insertToJmpMap();
            result.emplace_back(Sta{});
            break;

          case 5:
            //fprintf(f, "JMP\n");

            insertToJmpMap();
            result.emplace_back(Jmp{INT});
            break;

          case 6:
            //fprintf(f, "END\n");

            insertToJmpMap();
            result.emplace_back(End{});
            break;

          case 7:
            //fprintf(f, "RET\n");

            insertToJmpMap();
            result.emplace_back(Ret{});
            break;

          case 8:
            //fprintf(f, "DROP\n");

            insertToJmpMap();
            result.emplace_back(Drop{});
            break;

          case 9:
            //fprintf(f, "DUP\n");

            insertToJmpMap();
            result.emplace_back(Dup{});
            break;

          case 10:
            //fprintf(f, "SWAP\n");

            insertToJmpMap();
            result.emplace_back(Swap{});
            break;

          case 11:
            //fprintf(f, "ELEM\n");

            insertToJmpMap();
            result.emplace_back(Elem{});
            break;

          default:
            FAIL;
        }
        break;

      case 2:
      case 3:
      case 4: {
//        //fprintf (f, "%s\t\n", lds[h-2].c_str());
        insertToJmpMap();
        dsg::Designation d;
        switch (l) {
          case 0:
            //fprintf(f, "G()\n");
            d.emplace<dsg::Global>(dsg::Global{INT});
            break;
          case 1:
            //fprintf(f, "L()\n");
            d.emplace<dsg::Local>(dsg::Local{INT});
            break;
          case 2:
            //fprintf(f, "A()\n");
            d.emplace<dsg::Arg>(dsg::Arg{INT});
            break;
          case 3:
            //fprintf(f, "C()\n");
            d.emplace<dsg::Const>(dsg::Const{BOX(INT)});
            break;
          default:
            FAIL;
        }

        switch (h) {
          case 2:
            insertToJmpMap();
            result.emplace_back(Ld{d});
            break;
          case 3:
            insertToJmpMap();
            result.emplace_back((Lda{d}));
            break;
          case 4:
            insertToJmpMap();
            result.emplace_back((St{d}));
            break;
          default:
            FAIL;
        }
        break;
      }
      case 5:
        switch (l) {
          case 0:
            //fprintf(f, "CJMPz\n");

            insertToJmpMap();
            result.emplace_back((CJmp{false, INT}));
            break;

          case 1:
            //fprintf(f, "CJMPnz\n");

            insertToJmpMap();
            result.emplace_back((CJmp{true, INT}));
            break;

          case 2:
            //fprintf(f, "BEGIN\n");

            insertToJmpMap();
            result.emplace_back((Begin{INT, INT}));
            break;

          case 3:
            //fprintf(f, "CBEGIN\n");

            insertToJmpMap();
            result.emplace_back((CBegin{INT, INT}));
            break;

          case 4:
            //fprintf(f, "CLOSURE\n");
            /*insertToJmpMap();            *//*//fprintf(f, "CLOSURE\t0x%.8x", INT);
            {
              int n = INT;
              for (int i = 0; i < n; i++) {
                switch (BYTE) {
                  case 0:
                    //fprintf(f, "G(%d)", INT);
                    break;
                  case 1:
                    //fprintf(f, "L(%d)", INT);
                    break;
                  case 2:
                    //fprintf(f, "A(%d)", INT);
                    break;
                  case 3:
                    //fprintf(f, "C(%d)", INT);
                    break;
                  default:
                    FAIL;
                }
              }
            };*//*
            result.emplace_back((Closure{INT}));*/
            break;

          case 5:
            //fprintf(f, "CALLC\n");

            insertToJmpMap();
            result.emplace_back((CallC{INT}));
            break;

          case 6:
            //fprintf(f, "CALL\n");

            insertToJmpMap();
            result.emplace_back((Call{INT, INT}));
            break;

          case 7:
            //fprintf(f, "TAG\n");

            insertToJmpMap();
            result.emplace_back((Tag{STRING, INT}));
            break;

          case 8:
            //fprintf(f, "ARRAY\n");

            insertToJmpMap();
            result.emplace_back((Array{INT}));
            break;

          case 9:
            //fprintf(f, "FAIL\n");

            insertToJmpMap();
            result.emplace_back((Fail{INT, INT}));
            break;

          case 10:
            //fprintf(f, "LINE\n");

            insertToJmpMap();
            result.emplace_back((Line{INT}));
            break;

          default:
            FAIL;
        }
        break;

      case 6:
//        //fprintf (f, "PATT\t%s\n", pats[l].c_str());
//        insertToJmpMap();
        result.emplace_back((Patt{pats[l]}));
        break;

      case 7: {
        switch (l) {
          case 0:
            //fprintf(f, "CALL\tLread\n");

            insertToJmpMap();
            result.emplace_back(call::CallRead{});
            break;

          case 1:
            //fprintf(f, "CALL\tLwrite\n");

            insertToJmpMap();
            result.emplace_back(call::CallWrite{});
            break;

          case 2:
            //fprintf(f, "CALL\tLlength\n");

            insertToJmpMap();
            result.emplace_back(call::CallLength{});
            break;

          case 3:
            //fprintf(f, "CALL\tLstring\n");

            insertToJmpMap();
            result.emplace_back(call::CallString{});
            break;

          case 4:
            //fprintf(f, "CALL\tBarray\n");

            insertToJmpMap();
            result.emplace_back(call::CallArray{INT});
            break;

          default:
            FAIL;
        }
      }
        break;

      default:
        FAIL;
    }
  } while (true);

  stop:
  for (auto &token: result) {
    std::visit([&](auto &&arg) -> void {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<CJmp, T> || std::is_same_v<Jmp, T> || std::is_same_v<Call, T>) {
          arg.i = jmpMap.find(arg.i)->second;
        }
    }, token);
  }

  return result;
}

namespace stack {
    const int SIZE = 1024 * 1024;

    int array[SIZE];
    int top = SIZE;

    void push(int x) {
      if (__gc_stack_top <= (size_t) &stack::array) {
        throw std::overflow_error("Stack Overflow");
      }
      array[--top] = x;
      __gc_stack_top = (size_t) &array[top];
    }

    int pop() {
      if (__gc_stack_top == __gc_stack_bottom) {
        throw std::out_of_range("Stack Underflow");
      }
      int i = array[top++];
      __gc_stack_top = (size_t) &array[top];
      return i;
    }

    int peek() {
      if (__gc_stack_top == __gc_stack_bottom) {
        throw std::out_of_range("Stack Underflow");
      }
      return array[top];
    }
}

int main(int argc, char *argv[]) {
  bytefile *f = read_file(argv[1]);

  std::vector<Token> tokens = disassemble(f);

  __gc_stack_bottom = (size_t) &stack::array + sizeof(stack::array);
  __gc_stack_top = (size_t) &stack::array + sizeof(stack::array);
  __init();

  std::vector<Scope> scopes;
  std::vector<int> global_area(f->global_area_size);

  std::vector<int> rets;

#define IS(A, Type) constexpr (std::is_same_v<A, Type>)

  for (int i = 0; i < tokens.size(); ++i) {
    auto &token = tokens[i];
    std::visit([&](auto &&arg) -> void {
        using T = std::decay_t<decltype(arg)>;
        if IS(BinOp, T) {
          const auto y = UNBOX(stack::pop());
          const auto x = UNBOX(stack::pop());

          switch (arg.op) {
            case Operator::PLUS:
              stack::push(BOX(x + y));
              break;

            case Operator::MINUS:
              stack::push(BOX(x - y));
              break;

            case Operator::MULT:
              stack::push(BOX(x * y));
              break;

            case Operator::DIV:
              stack::push(BOX(x / y));
              break;

            case Operator::MOD:
              stack::push(BOX(x % y));
              break;

            case Operator::LT:
              stack::push(BOX(x < y));
              break;

            case Operator::LE:
              stack::push(BOX(x <= y));
              break;

            case Operator::GT:
              stack::push(BOX(x > y));
              break;

            case Operator::GE:
              stack::push(BOX(x >= y));
              break;

            case Operator::EQ:
              stack::push(BOX(x == y));
              break;

            case Operator::NEQ:
              stack::push(BOX(x != y));
              break;

            case Operator::DISJ:
              stack::push(BOX(x || y));
              break;

            case Operator::CONJ:
              stack::push(BOX(x && y));
              break;

            default:
              break;
          }
        } else if IS(Const, T) {
          stack::push(arg.x);
        } else if IS(Sexp, T) {
          auto tag = LtagHash(arg.s);
          auto bsexpTag = BsexpTag(BOX(arg.n + 1), tag);
          for (int j = 0; j < arg.n; ++j) {
            Bsta((void *) stack::pop(), BOX(arg.n - j - 1), bsexpTag);
          }

          stack::push(int(bsexpTag));
        } else if IS(Sti, T) {
          auto value = stack::pop();
          auto ref = stack::pop();
          *(int *) ref = value;
          stack::push(value);
        } else if IS(Sta, T) {
          auto v = (void *) stack::pop();
          auto ind = stack::pop();
          auto x = (void *) stack::pop();

          stack::push(int(Bsta(v, ind, x)));
        } else if IS(Jmp, T) {
          i = arg.i - 1;
        } else if IS(CJmp, T) {
          auto top = UNBOX(stack::pop());

          if (arg.notZero && top != 0 || !arg.notZero && top == 0) {
            i = arg.i - 1;
          }
        } else if IS(Begin, T) {
          auto scope = Scope{std::vector<int>(arg.n), std::vector<int>(arg.m)};
          scopes.push_back(std::move(scope));

          if (!rets.empty()) {
            for (size_t j = 0; j < arg.n; ++j) {
              scopes.back().args[arg.n - j - 1] = stack::pop();
            }
          }
        } else if IS(End, T) {
          scopes.pop_back();
          if (rets.empty()) {
            i = tokens.size() - 1;
          } else {
            i = rets.back();
            rets.pop_back();
          }
        } else if IS(Call, T) {
          rets.push_back(i);
          i = arg.i - 1;
        } else if IS(Drop, T) {
          stack::pop();
        } else if IS(Dup, T) {
          stack::push(stack::peek());
        } else if IS(Swap, T) {
          auto a = stack::pop();
          auto b = stack::pop();

          stack::push(a);
          stack::push(b);
        } else if IS(Elem, T) {
          auto ind = stack::pop();
          auto p = (void *) stack::pop();

          stack::push(int(Belem(p, ind)));
        } else if IS(Ld, T) {
          std::visit([&](auto &&d) -> void {
              using D = std::decay_t<decltype(d)>;
              if IS(dsg::Const, D) {
                stack::push(d.x);
              } else if IS(dsg::Global, D) {
                stack::push(global_area[d.x]);
              } else if IS(dsg::Local, D) {
                stack::push(scopes.back().locals[d.x]);
              } else if IS(dsg::Arg, D) {
                stack::push(scopes.back().args[d.x]);
              }
          }, arg.d);
        } else if IS(Lda, T) {
          std::visit([&](auto &&d) -> void {
              using D = std::decay_t<decltype(d)>;
              if IS(dsg::Global, D) {
                stack::push(global_area[d.x]);
              } else if IS(dsg::Local, D) {
                stack::push(scopes.back().locals[d.x]);
              } else if IS(dsg::Arg, D) {
                stack::push(scopes.back().args[d.x]);
              }
          }, arg.d);
        } else if IS(St, T) {
          std::visit([&](auto &&d) -> void {
              using D = std::decay_t<decltype(d)>;
              if IS(dsg::Global, D) {
                global_area[d.x] = stack::peek();
              } else if IS(dsg::Local, D) {
                scopes.back().locals[d.x] = stack::peek();
              } else if IS(dsg::Arg, D) {
                scopes.back().args[d.x] = stack::peek();
              }
          }, arg.d);
        } else if IS(Tag, T) {
          stack::push(Btag((void *) stack::pop(), LtagHash(arg.s), BOX(arg.n)));
        } else if IS(Array, T) {
          stack::push(Barray_patt((void *) stack::pop(), BOX(arg.n)));
        } else if IS(call::Call, T) {
          std::visit([&](auto &&c) -> void {
              using C = std::decay_t<decltype(c)>;
              if IS(call::CallRead, C) {
                stack::push(Lread());
              } else if IS(call::CallWrite, C) {
                Lwrite(stack::pop());
                stack::push(BOX(0));
              } else if IS(call::CallArray, C) {
                auto x = LmakeArray(BOX(c.n));
                for (size_t i = 0; i < c.n; ++i) {
                  Bsta((void *) stack::pop(), BOX(c.n - i - 1), x);
                }
                stack::push(int(x));
              } else if IS(call::CallLength, C) {
                auto res = Llength((void *) stack::pop());
                stack::push(res);
              }
          }, arg);
        }
    }, token);
  }

  return 0;
}