# Contributing

Sapphire++ aims to become a community project with use cases and members from a
wide range of research fields. It is our goal to build an inclusive and
participatory community, so we are happy that you are interested in
participating. To this end, we have adopted a
[code of conduct](https://sapphirepp.org/latest/md_doc_2pages_2code-of-conduct.html)
that we expect all participants to adhere to.


## Getting started with git and GitHub

If you are new to using git or the GitHub platform you may find it helpful to
review some resources available on the web.

@todo Add our favourite git and GitHub intros here. Maybe also link
 [lecture 32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html).


## Getting help

There is a Sapphire++ mailing list to get in contact with the developers and our
community. It is the perfect place for questions and discussions about the
library on all levels. Feel free to
[subscribe](https://mein.manitu.de/public/webhosting/mailinglist/?id=156396&auth=tSGpYMy4VrSEjX9vtxJFpMsgFDbfjT1a).

The mailing list is a public forum, which can be accessed under
[https://mailinglist.sapphirepp.org](https://mailinglist.sapphirepp.org) and
which we consider a resource for any kind of problems users encounter when they
use Sapphire++.

We note that all emails sent to the list are distributed both to the list
subscribers and copied to this public archive, for people to browse or search
without the need to be subscribed.

Subscribing to the mailing list is equivalent to confirming that you are okay
with publishing the emails you sent to the mailing list on this archive.
Nonetheless, we use all possibilities, which the mailing list archive software
provides, to hide sensible and irrelevant information. For example, your email
address will be disguised.

Please do not post any confidential information (e.g. email addresses of others)
and follow our
[code of conduct](https://sapphirepp.org/latest/md_doc_2pages_2code-of-conduct.html).

To unsubscribe follow this
[link](https://mein.manitu.de/public/webhosting/mailinglist/?id=156396&auth=tSGpYMy4VrSEjX9vtxJFpMsgFDbfjT1a)
and enter your email address. If you would like us to delete any of your
contributions, please write an email to
[contact@sapphirepp.org](mailto:contact@sapphirepp.org).


## Bug reports

It is a great help to us if you report any bugs in the library that you may
find. We keep track of all open issues
[here](https://github.com/sapphirepp/sapphirepp/issues). If you can, please try
to include a minimal failing example that can help us reproduce the problem.


## Making changes and contributing new features

To make a change to Sapphire++ you should create a *fork* (through GitHub) and a
separate *branch* (sometimes called a feature branch). You can propose that your
branch be combined with the main branch by opening a *pull request*. This will
give a chance for others to review your code. While this seems very formal,
keeping the code review in one place makes it easier to coordinate changes.
Please do not hesitate to ask questions about the workflow on the mailing list
if you are not sure what to do.

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
