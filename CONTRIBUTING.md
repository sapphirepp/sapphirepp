# Contributing to Sapphire++

Sapphire++ aims to become a community project with use cases and members from a
wide range of research fields. It is our goal to build an inclusive and
participatory community, so we are happy that you are interested in
participating. To this end, we have adopted a [code of
conduct](https://sapphirepp.org/latest/md_doc_2pages_2code-of-conduct-page.html)
that we expect all participants to adhere to.×

## Getting started with git and GitHub

If you are new to using git or the GitHub platform you may find it helpful to
review some resources available on the web. For example, the [Pro Git
Book](https://git-scm.com/book/en/v2) by Scott Chacon and Ben Straub or the
[deal.II video lectures](https://www.math.colostate.edu/~bangerth/videos.html)
[32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html) and
[32.75](https://www.math.colostate.edu/~bangerth/videos.676.32.75.html) by
Wolfgang Bangerth.

## Getting help

We use
[GitHub Discussions](https://github.com/sapphirepp/sapphirepp/discussions)
as a forum for discussions.
It is the perfect place for questions
and discussions about the library on all levels.
It also serves as a resource for any kind of problems users encounter
when they use Sapphire++.

To ask a question or comment,
you need to have a
[GitHub Account](https://docs.github.com/en/get-started/signing-up-for-github/signing-up-for-a-new-github-account).
To stay updated,
you can either subscribe to individual discussions,
or [watch the repository](https://docs.github.com/en/github/managing-subscriptions-and-notifications-on-github/setting-up-notifications/about-notifications)
to receive notifications about all new announcements and developments.

## Bug reports

It is a great help to us if you report any bugs in the library that you may
find. We keep track of all open issues
[here](https://github.com/sapphirepp/sapphirepp/issues). If you can, please try
to include a minimal failing example that can help us reproduce the problem.

## Making changes and contributing new features

To make a change to Sapphire++ you should create a *fork* (through GitHub)
and a separate *branch* (sometimes called a feature branch).
You can propose that your branch be combined with the main branch
by opening a *pull request*.
This will give a chance for others to review your code.
While this seems very formal,
keeping the code review in one place makes it easier to coordinate changes.
Please do not hesitate to ask questions about the workflow
on the GitHub Discussions page if you are not sure what to do.

Since Sapphire++ builds heavily on the [deal.II](https://www.dealii.org) library
for finite elements, we recommand to familiarize yourself with the
[deal.II tutorials](https://www.dealii.org/current/doxygen/deal.II/Tutorial.html)
or the
[deal.II lectures](https://www.math.colostate.edu/~bangerth/videos.html).
However, Sapphire++ is build such that you can use it without knowledge of
[deal.II](https://www.dealii.org). If you want to contribute new use cases or
examples on a top level, you may skip this step.

To keep the style of the source code consistent we use a set of
[coding conventions](https://sapphirepp.org/latest/coding-conventions.html).
This convention essentially follows the
[deal.II coding conventions](https://www.dealii.org/developer/doxygen/deal.II/CodingConventions.html)
by using `clang-format` for indentation, camel case for classes, and lower case
letters with underscores for everything else. If you are new to the project then
we will work with you to ensure your contributions are formatted with this
style, so please do not think of it as a road block if you would like to
contribute some code.

Sapphire++ is licensed under the GNU Lesser General Public License; while you
will retain copyright on your contributions, all changes to the core library
must be provided under this common license.
