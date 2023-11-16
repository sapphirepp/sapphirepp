# Coding conventions {#coding-conventions}

@tableofcontents

Throughout Sapphire++, we strive to keep our programming style and the kind of
interfaces we provide as consistent as possible. To this end, we have adopted a
set of coding conventions that we attempt to follow wherever possible. They have
two parts: style issues, and something we call "defensive programming", the
latter being an attempt to let our code help us find bugs. When reading through
them, it is important to remember that styles are not god given or better than
any other set of conventions; their purpose is merely to keep Sapphire++ as
uniform as possible. Uniformity reduces the number of bugs we produce because we
can, for example, always assume that input arguments come before output
arguments of a function call. They also simplify reading code because some
things become clear already by looking at the style a piece of code is written,
without having to look up the exact definition of something.

These coding conventions are based on the
[deal.II coding conventions](https://www.dealii.org/developer/doxygen/deal.II/CodingConventions.html)
with only a few modifications.

@note These coding conventions are to be understood as goal. You may find, that
  we not always follows these conventions in the code yet. You are welcome to
  help us to improve our code base by contributing style improvements!


## Notes indentation {#coding-conventions-indentation}

We use clang-format to format our code. We recommand to integrate code
formatting in your IDE.


## Style issues {#coding-conventions-style}

1. Traditional logical operators should be used instead of their English
   equivalents (i.e., use &&, ||, and ! instead of and, or, and not).

2. In the implementation files, after each function, three empty lines are
   expected to enable better readability. One empty line occurs in functions to
   group blocks of code, since two empty lines are not enough to visibly
   distinguish sufficiently that the code belongs to two different functions.

3. Whenever an integer variable can only assume non-negative values, it is
   marked as unsigned. The same applies to functions that can only return
   positive or zero values. Example: `unsigned int dim`.

4. Whenever an argument to a function will not be changed, it should be marked
   `const`, even if passed by value. Generally, we mark input parameters as
   `const`. This aids as an additional documentation tool to clarify the intent
   of a parameter (input, output, or both) and lets the compiler issue warnings
   if such a parameter is changed, which is often either involuntarily or poor
   style.

5. Whenever a function does not change any of the member variable of the
   embedding class/object, it should be marked as `const`.

6. Function and variable names may not consist of only one or two letters,
   unless the variable is a pure counting index.

7. Type aliases (`using`-declarations) are preferred to `typedef-declarations`.

8. The layout of class declarations is the following: first the block of public
   functions, beginning with the constructors, then the destructors. If there
   are public member variables, these have to occur before the constructor.
   Public variables shall only be used if constant (in particular if they are
   static and constant) or unavoidable (*Exception: `Parameters` classes*).  
   After the public members, the protected and finally the private members are
   to be listed. The order is as above: first variables then functions.  
   Exceptions shall be declared at the end of the public section before the
   non-public sections start.  
   We do not use the C++11-style class member initialization for member
   variables that are neither `static const` nor `static constexpr`; i.e.,
   instead of  

    ```cpp
    class Foo
    {
      int a = 42;
      int *b = nullptr;
    };
    ```

    write

    ```cpp
    class Foo
    {
      Foo();

      int a;
      int *b;
    };
    

    inline Foo::Foo()
    : a(42)
    , b(nullptr)
    {}
    ```

9. Runtime parameters are always to be handled in `*Parameters` classes. The run
   time variables are to be declared as **non-constant** `public` member variables
   (non-constant because they are set at run time and `public` for easy access).
   They are declared in `declare_parameters()` and parsed in
   `parse_parameters()`. All other functions and classes should only ever get
   `const` access to the `Parameters` classes. Example:
   @vfpref{VFPParameters}.

10. If a function has both input and output parameters, usually the input
    parameters shall precede the output parameters, unless there are good
    reasons to change this order. (The most common reason is trailing input
    parameters with default values.)

11. Exceptions are used for internal parameter checking and for consistency
    checks through the @dealii `Assert` macro. Exception handling like done by
    the C++ language (`try/throw/catch`, and using the `AssertThrow` macro) are
    used to handle run time errors (like I/O failures) which must be on in any
    case, not only in debug mode.

12. Sometimes it makes sense to implement a class by using several non-member
    functions that are not part of the public interface and are only meant to be
    called in the current source file. Such free functions should be put in an
    internal namespace structured in the following way:

    ```cpp
    namespace internal
    {
      namespace ClassNameImplementation
      {
        // free functions go here
      }
    }
    ```

    where `ClassName` is the name of the calling class.

13. Classes, namespaces and types generally are named using uppercase letters to
    denote word beginnings (e.g. `OutputParameters`) — sometimes called
    *[camel case](https://en.wikipedia.org/wiki/Camel_case)* — while functions
    and variables use lowercase letters and underscores to separate words.

14. We may use forward declarations in header files to, hopefully, improve
    compilation speeds by not using headers when we just need to mark a certain
    type as an argument to a function. The convention is that, if all we need is
    a type name, then the type may be forward declared in the header where we
    need it; if a function (or member function) can return a value then a
    declaration of that value's type should be available (by including the
    necessary header).


## Instantiation of templated functions/classes {#coding-conventions-templates}

The majority of classes and functions in @sapphire are templated. This brings a
question of how and where such objects are instantiated, if at all. We adopt the
following convention:

1. If we can enumerate all possible template arguments (e.g., the dimension can
   only be 1, 2, or 3), then a function template goes into the .cc file, and we
   explicitly instantiate all possibilities. Users will not have any need to
   ever see these function templates because they will not want to instantiate
   these functions for any other template arguments anyway.

2. If we can not enumerate all possible template arguments (e.g., vector types –
   because users might want to define their own vector kinds) but at least know
   a few common usage cases, then the function is put into a `.templates.h`
   file. We `#include` it into the `.cc` file and instantiate the functions for
   all the common arguments. For almost all users, this will be just fine – they
   only use the (vector, matrix, ...) types we already instantiate, and for them
   the `.templates.h` file will not be of any interest. It will also not slow
   down their compilations because nothing they see will `#include` the
   `.templates.h` file. But users who define their own (vector, matrix, ...)
   types can instantiate the template functions with their own user-defined
   types by including the `.templates.h` files.

3. Finally, if we can not assume in advance which values template arguments will
   take, the definitions of functions are provided at the bottom of the header
   file with declarations. The definitions should be guarded with
   `#ifndef DOXYGEN ... #endif` to prevent Doxygen from picking them up.


## Defensive programming {#coding-conventions-defensive-programming}

Defensive programming is a term that we use frequently when we talk about
writing code while in the mindset that errors will happen. Here, errors can come
in two ways: first, I can make a mistake myself while writing a function; and
secondly, someone else can make a mistake while calling my function. In either
case, I would like to write my code in such a way that errors are (i) as
unlikely as possible, (ii) that the compiler can already find some mistakes, and
(iii) that the remaining mistakes are relatively easy to find, for example
because the program aborts. Defensive programming is then a set of strategies
that make these goals more likely.

Over time, we have learned a number of techniques to this end, some of which we
list here:

1. *Assert preconditions on parameters:* People call functions with wrong or
   nonsensical parameters, all the time. As the prototypical example, consider a
   trivial implementation of vector addition:

    ```cpp
    Vector &
    operator+=(Vector       &lhs,
               const Vector &rhs)
    {
      for (unsigned int i=0; i<lhs.size(); ++i)
        lhs(i) += rhs(i);
      return lhs;
    }
    ```

    While correct, this function will get into trouble if the two vectors do not
    have the same size. You think it is silly to call this function with vectors
    of different size? Yes, of course it is. But it happens all the time: people
    forget to reinitialize a vector, or it is reset in a different function,
    etc. It happens. So if you are in such an unlucky case, it can take a long
    time to figure out what's going on because you are likely to just read
    uninitialized memory, or maybe you are writing to memory the `lhs` vector
    doesn't actually own. Neither is going to lead to immediate termination of
    the program, but you'll probably get random errors at a later time. It would
    be much easier if the program just stopped here right away. The following
    implementation will do exactly this:

    ```cpp
    Vector &
    operator+=(Vector       &lhs,
               const Vector &rhs)
    {
      AssertDimension (lhs.size(), rhs.size());
      for (unsigned int i=0; i<lhs.size(); ++i)
        lhs(i) += rhs(i);
      return lhs;
    }
    ```

    The `Assert` macro ensures that the condition is true at run time, and
    otherwise prints a string containing information encoded by the second
    argument and aborts the program. This way, when you write a new program that
    happens to call this function, you will learn of your error right away and
    have the opportunity to fix it without ever having to seriously debug
    anything.

    As a general guideline, whenever you implement a new function, think about
    the *preconditions* on parameter, i.e. what does the function expect to be
    true about each of them, or their combination. Then write assertions for all
    of these preconditions. This may be half a dozen assertions in some cases
    but remember that each assertion is a potential bug already found through
    trivial means.

    In a final note, let us remark that assertions are of course expensive: they
    may make a program 3 or 5 times slower when you link it against the debug
    version of the library. But if you consider your *overall* development time,
    the ability to find bugs quickly probably far outweighs the time you spend
    waiting for your program to finish. Furthermore, calls to the `Assert` macro
    are removed from the program in optimized mode (which you presumably only
    use once you know that everything runs just fine in debug mode). The
    optimized libraries are faster by a factor of 3-5 than the debug libraries,
    at the price that it's much harder to find bugs.  
    We highly recommend, to use the debug mode to develop and test your code,
    but to use the optimized mode for simulation runs.

2. *Assert postconditions:* If a function computes something non-trivial there
   may be a bug in the code. To find these, use postconditions: just like you
   have certain knowledge about useful values for input parameters, you have
   knowledge about what you expect possible return values to be. For example, a
   function that computes the norm of a vector would expect the norm to be
   positive. You can write this as follows:

    ```cpp
    double norm(const Vector &v)
    {
      double s = 0;
      for (unsigned int i=0; i<v.size(); ++i)
        s += v(i) * v(i);

      Assert (s >= 0, ExcInternalError());
      return std::sqrt(s);
    }
    ```

    This function is too simple to really justify this assertion, but imagine
    the computation to be lengthier, and you can see how the assertion helps you
    ensure (or hedge) yourself against mistakes. Note that one could argue that
    the assertion should be removed once we've run the program a number of times
    and found that the condition never triggers. But it's better to leave it
    right where it is: it encodes for the future (and for readers) knowledge you
    have about the function; if someone comes along and replaced the
    implementation of the function by a more efficient algorithm, the assertion
    can help make sure that the function continues to do what it is supposed to
    do.

3. *Assert internal states:* In a similar vein, if you have a complex algorithm,
   use assertions to ensure that your mental model of what is going on matches
   what is indeed true. For example, assume you are writing a function that
   ensures that mesh sizes do not change too much locally. You may end up with a
   code of the following kind:

    ```cpp
    for (const auto &cell = triangulation.active_cell_iterators())
      for (unsigned int face=0; ...)
        {
          if (something)
            { ... }
          else
            {
              // we have a cell whose neighbor must
              // be at the boundary if we got here
            }
        }
    ```

    The conditions that got us into the else-branch may be complicated, and
    while it may be true that we believed that the only possibility we got here
    is that the neighbor is at the boundary, there may have been a bug in our
    implementation. There may also have been a bug in our thinking, or someone
    changes the code way above in the same function and forgets about the issue
    here, or a change at a completely different location in the library makes
    the assumption untenable. In all of these cases, the explicit statement of
    our assertion makes sure that these problems are easily found.

4. *Initialize variables at the point of their declaration if they live on the
   stack:* Traditional C required that variables are declared at the beginning
   of the function even if they are only used further below. This leads to code
   like this that we may imagine in a 1d code:

    ```cpp
    template <unsigned int dim>
    void foo ()
    {
      Point<dim> cell_center;
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
          ...
        }
      ...
    }
    ```

    The problem is that if the code between the declaration and initialization
    is long and complicated, you can't look up on one page what the type of a
    variable is and what it's value may be. In fact, it may not even be quite
    clear that the variable is used initialized at all, or whether it is
    accidentally left uninitialized.

    A better way to do this would be as follows:

    ```cpp
    template <unsigned int dim>
    void foo ()
    {
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          Point<dim> cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
          ...
        }
      ...
    }
    ```

    This makes it much clearer what the type of the variable is and that it is
    in fact only ever used when initialized. Furthermore, if someone wants to
    read the code to see what the variable is in fact doing, declaring and
    initializing it in the innermost possible scope makes this task easier: we
    don't have to look upwards for it beyond the declaration, and we don't have
    to look downward beyond the end of the current scope since this is where the
    variable dies.

    As a final note, it is clear that you can only do this sort of stuff for
    variables that completely live on the stack without allocating memory on the
    heap. This is only true for builtin types like `int`, `double`, `char`, etc,
    as well as the @dealref{Point} and @dealref{Tensor} classes. Everything else
    has something like a `std::vector` as a member variable, which requires
    memory allocation — you don't want to declare these inside loops, at least
    not if the loop is traversed frequently.

5. *Make variables `const`:* To pick up on the example above, note that in most
   cases we will never change the variable so initialized any more. In other
   words, if this is the case, we may as well write things as follows:

    ```cpp
    template <unsigned int dim>
    void foo ()
    {
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          const Point<dim> cell_center = (cell->vertex(0) +
                                          cell->vertex(1)) / 2;
          ...
        }
      ...
    }
    ```

    By marking the variable as constant we make sure that we don't accidentally
    change it. For example, the compiler could catch code like this:

    ```cpp
    if (cell_center[0] = 0)
      ...
    ```

    This was most likely meant to be a `==` rather than an assignment. By
    marking the variable as `const`, the compiler would have told us about this
    bug. Maybe equally importantly, human readers of the code need not look
    further down whether the value of the variable may actually be changed
    somewhere between declaration and use — it can't be if it is marked as
    `const`.

6. *Make input arguments of functions `const`:* The same essentially holds true
   as well as for function arguments: If you have no intention of changing a
   variable (which is typically the case for input arguments), then mark it as
   constant. For example, the following function should take its argument as a
   constant value:

   @todo Find better example, maybe `PDESystem::create_lms_matrix`?

    ```cpp
    template <unsigned int dim>
    typename Triangulation<dim>::cell_iterator
    CellAccessor<dim>::child(const unsigned int child_no)
    {
      ...
      return something;
    }
    ```

    Here, the user calls `cell->child(3)`, for example. There really is no
    reason why the function would ever want to change the value of the
    `child_no` argument — so mark it as constant: this both helps the reader of
    the code understand that this is an input argument of the function for which
    we need not search below whether it is ever changed, and it helps the
    compiler help us find bugs if we ever accidentally change the value.


## Documentation style {#coding-conventions-documentation}

We build the documentation for @sapphire using
[Doxygen](https://www.doxygen.nl/index.html). We use the `/**` style for comment
blocks, `///` for inline documentation and `@` for Doxygen commands. You can
[Markdown](https://www.doxygen.nl/manual/markdown.html) syntax in your comments
and include Latex formulas using `\f$ \alpha \f$` for inline and

```plain
\f[
  \alpha
\f]
```

for display style formulas.

To document classes and functions, use the `@brief` tag for a short summary of
the function or class, and the `@param` and `@return` tag to document each
parameter and the return value respectively. An example documentation of a
function would look like this:

```plain
  /**
   * @brief Add two numbers.
   *
   * This function add two numbers \f$ a \f$ and \f$ b \f$ and returns the
   * sum \f$ a + b \f$.
   *
   * @param a The first number \f$ a \f$.
   *        Indent following lines, if needed.
   * @param b The second number \f$ b \f$
   * @return double The sum \f$ a + b \f$
   */
  double
  add_numbers(const double a, const double b);
```

Every file starts with a header, first including the Copyright statement,
followed by a brief description of the file:

```plain
COPYRIGHT STATEMENT

/**
 * @file filename
 * @author Author 1 (email)
 * @author Author 2 (email)
 * @author ...
 * @brief A brief description of the contents of the file
 */

 #include ...
```

To acknowledge your contribution, we will add you to the author list of all
files you have contributed to. If you do not want to be listed, please let us
know!

High level documentation, like the introduction and examples, are written in
Markdown.  
We use some custom commands like `@ sapphire` and `@ dealii` (without the space)
to refer to @sapphire and @dealii consistently.


## Git {#coding-conventions-git}

We would like to keep a clean and linear git history. To this end, we would like
to ask you to `rebase` your branch on top of `main` or `devel` before opening a
pull request. We prefer `rebase` over `merge` because it keeps the history
linear. We also advice to use the following git configuration:

```bash
git config --local pull.rebase true
```

We use [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) to
ensure meaning full and consistent commit messages. Please follow these
conventions by using the following commit message format:

```text
type(scope): subject

optional body

optional footer
```

with one of the following type: `feat`, `fix`, `docs`, `style`, `refactor`,
`perf`, `test`, ...


Don't think of these conventions as a hurdle to contribute! Feel free to
contribute in any way you like. We will help you to get your contribution into
the right shape. We are happy to help you to get started!
