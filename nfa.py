#!/usr/bin/env python3

from nfa_regex import *
from state import *

try:
    from graphviz import Digraph
except:
    pass

class NFA:
    """A class for NFAs that can be constructed from a Regex"""

    def __init__(self, r=None):
        """Constructs an empty-language NFA when r==None, 
            or constructs an NFA from the provided Regex object."""
        self.r = r  # Save a reference to the regex itself
        if r is None:
            # If no Regex instance is provided: the empty language 
            start = State()
            self.Q = {start}
            self.Sigma = set()
            self.s = start
            self.F = set()
            self.delta = set()
        elif isinstance(r, EpsilonRegex):
            start = State()
            end = State()
            self.Q = {start, end}
            self.Sigma = set()
            self.s = start
            self.F = {end}
            self.delta = {
                (start, "", end)
            }
            pass  # DONE: the epsilon language, base form
        elif isinstance(r, StarRegex):
            nfa0 = NFA(r.r0)
            start = State()
            end = State()
            self.Q = {start, end} | nfa0.Q
            self.Sigma = nfa0.Sigma
            self.s = start
            self.F = {end}
            self.delta = {
                             (start, "", end),
                             (start, "", nfa0.s),
                             (min(nfa0.F), "", end),
                             (min(nfa0.F), "", nfa0.s)
                         } | nfa0.delta
            pass  # DONE: handle derived form star
        elif isinstance(r, SeqRegex):
            nfa0 = NFA(r.r0)
            nfa1 = NFA(r.r1)
            self.Q = nfa0.Q | nfa1.Q
            self.Sigma = nfa0.Sigma | nfa1.Sigma
            self.s = nfa0.s
            self.F = nfa1.F
            self.delta = {
                             (min(nfa0.F), "", nfa1.s)
                         } | nfa0.delta | nfa1.delta
            pass  # DONE: handle derived form juxtaposition/sequencing
        elif isinstance(r, DisjRegex):
            nfa0 = NFA(r.r0)
            nfa1 = NFA(r.r1)
            start = State()
            end = State()

            self.Q = {start, end} | nfa0.Q | nfa1.Q
            self.Sigma = nfa0.Sigma | nfa1.Sigma
            self.s = start
            self.F = {end}
            self.delta = ({(start, "", nfa0.s), (start, "", nfa1.s),
                           (min(nfa0.F), "", end), (min(nfa1.F), "", end)}
                          | nfa0.delta | nfa1.delta)
        elif isinstance(r, CharRegex):
            start = State()
            end = State()
            self.Q = {start, end}
            self.Sigma = {str(r)}
            self.s = start
            self.F = {end}
            self.delta = {
                (start, str(r), end)
            }
            pass  # DONE: handle base form
        else:
            raise "NFA must be constructed from a Regex (or None)."

    def matches(self, s):
        # simulate an NFA or DFA (this code can work for both!)
        current = self.epsilonClosure({self.s})

        for x in s:
            current = self.epsilonClosure(self.move(current, x))
        if current & self.F == set():
            return False
        else:
            return True

    def dfs(self, qs, ext, c):
        ext.add(qs)
        for state in self.Q:
            if state not in ext and (qs, c, state) in self.delta:
                self.dfs(state, ext, c)

    def epsilonClosure(self, qs):
        """Returns the set of states reachable from those in qs without consuming any characters."""
        # pass  # DONE: write epsilon closure (see slides)
        res = set()
        for q in qs:
            self.dfs(q, res, "")
        return res

    def move(self, qs, x):
        """Returns the set of states reachable from those in qs by consuming character x."""
        # pass  # DONE: write move (see slides)
        res = set()
        for q in qs:
            for state in self.Q:
                if (q, x, state) in self.delta:
                    self.dfs(state, res, x)
        return res

    def NFA_to_DFA(self):
        """Returns a DFA equivalent to self."""
        dfa = NFA()
        dfa.r = self.r
        start_set = self.epsilonClosure({self.s})
        start = State(frozenset(start_set))
        dfa.Q = {start}
        # DONE: Construct dfa from self (using subset algorithm)
        dfa.Sigma = self.Sigma
        dfa.s = start
        dfa.delta = set()
        ext = set()
        diff = dfa.Q.difference(ext)
        while diff:
            ext |= diff
            for c in dfa.Sigma:
                for dstate in diff:
                    current = self.epsilonClosure(self.move(dstate.name, c))
                    current_state = State(frozenset(current))
                    if current_state not in dfa.Q:
                        dfa.Q.add(current_state)
                    dfa.delta |= {(dstate, c, current_state)}
            diff = dfa.Q.difference(ext)
        dfa.F = set()
        for dstate in dfa.Q:
            if dstate.name & self.F:
                dfa.F.add(dstate)
        return dfa

    def statecount(self):
        return len(self.Q)

    def isDFA(self):
        # Checks to see if the NFA is also a DFA:
        # i.e., reports True iff self.delta is a function
        for q in self.Q:
            outgoingset = set()
            for e in self.delta:
                if q == e[0]:
                    if e[1] in outgoingset or e[1] == "":
                        return False
                    outgoingset.add(e[1])
        return True

    def minimizeDFA(self):
        """Takes a DFA and returns a minimized DFA"""
        if not self.isDFA():
            raise "minimizeDFA must be provided a DFA"

        # Use partition refinement (to a fixpoint)
        parts = set()
        work = set()
        if self.F != set():
            parts.add(frozenset(self.F))
            work.add(frozenset(self.F))

        parts.add(frozenset(self.Q - self.F))
        mdfa = NFA()
        mdfa.Sigma = self.Sigma
        mdfa.r = self.r
        while work:
            p = work.pop()
            for c in self.Sigma:
                preds = set()
                for t in self.delta:
                    if t[1] == c and t[2] in p:
                        preds.add(t[0])

                for part in parts:
                    intersect = part & preds
                    diff = part - preds
                    if intersect and diff:
                        parts = (parts - {part}) | {intersect, diff}

                        if part in work:
                            work = (work - {part}) | {intersect, diff}
                        else:
                            if len(intersect) <= len(diff):
                                work.add(intersect)
                            else:
                                work.add(diff)
                        break

        mdfa.Q = set([State(x) for x in parts])
        mdfa.F = set([State(x) for x in parts if (self.F & x)])
        mdfa.s = list(filter(lambda x: self.s in x.name, mdfa.Q))[0]
        for t in self.delta:
            for start_t in mdfa.Q:
                if t[0] in start_t.name:
                    for end_t in mdfa.Q:
                        if t[2] in end_t.name:
                            mdfa.delta.add((start_t, t[1], end_t))
        return mdfa

    # generating visualizations of the NFA/DFA for understanding and debugging purposes
    def generateSVG(self, file_name="nfa", title=True):
        """Writes the current NFA to a dot file and runs graphviz (must be locally installed!)"""

        # Setup
        dot = Digraph(name='nfa', comment='nfa')
        dot.attr(rankdir='LR')
        names = {}
        if title:
            dot.node("regex", label='<<FONT POINT-SIZE="24">' + str(self.r) + '</FONT>>', shape="square",
                     style="rounded", height="0", width="0", margin="0.05")
        elif isinstance(title, str):
            dot.node("regex", label='<<FONT POINT-SIZE="24">' + title + '</FONT>>', shape="square", style="rounded",
                     height="0", width="0", margin="0.05")
        dot.node("*", style="invis", height="0", width="0", margin="0")

        def pad_lab(s):
            """For some reason graphviz needs spaces padding this to render right."""
            return "  " + s + "  "

        def str_lab(s):
            if s == "":
                return "ε"
            else:
                return str(s)

        # Nodes and Edges
        def namer(n):
            return "<q<sub><font point-size=\"11\">" + n + "</font></sub>>"

        if self.s in self.F:
            dot.node(str(len(names)), namer(str(len(names))), shape="doublecircle")
        else:
            dot.node(str(len(names)), namer(str(len(names))), shape="circle")
        names[self.s] = len(names)
        for n in (self.Q - self.F - {self.s}):
            dot.node(str(len(names)), namer(str(len(names))), shape="circle")
            names[n] = len(names)
        for n in (self.F - {self.s}):
            dot.node(str(len(names)), namer(str(len(names))), shape="doublecircle")
            names[n] = len(names)
        pseudodelta = dict()
        for e in self.delta:
            if (e[0], e[2]) in pseudodelta:
                pseudodelta[(e[0], e[2])] |= frozenset({e[1]})
            else:
                pseudodelta[(e[0], e[2])] = frozenset({e[1]})
        for k in pseudodelta.keys():
            dot.edge(str(names[k[0]]), str(names[k[1]]), label=pad_lab(
                ",".join(list(map(str_lab, sorted(list(pseudodelta[k])))))))
        dot.edge("*", str(names[self.s]), label="")

        dot.format = 'svg'
        dot.render(filename=file_name, cleanup=True)
