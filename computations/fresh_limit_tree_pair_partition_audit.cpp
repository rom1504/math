// Exact leading-diagram audit. Integer arithmetic only; no signing is sampled.
// Build: c++ -O3 -std=c++17 this_file.cpp -o /home/math/quadra/tmp/tree_pair_audit
// Run: /home/math/quadra/tmp/tree_pair_audit

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using U64 = std::uint64_t;

struct Tree {
    std::string name, canonical;
    std::vector<int> parent;
    std::vector<int> children;
    U64 automorphisms = 1;
};

std::vector<Tree> trees;

int add_tree(const std::string& name, std::vector<int> children) {
    Tree t;
    t.name = name;
    t.children = children;
    t.parent.push_back(-1);
    std::vector<std::string> child_names;
    std::map<std::string,int> multiplicity;
    for (int child: children) {
        const Tree& c = trees[child];
        int offset = static_cast<int>(t.parent.size());
        for (int p: c.parent) t.parent.push_back(p < 0 ? 0 : p+offset);
        child_names.push_back(c.canonical);
        multiplicity[c.canonical]++;
        t.automorphisms *= c.automorphisms;
    }
    std::sort(child_names.begin(),child_names.end());
    t.canonical = "(";
    for (const auto& s: child_names) t.canonical += s;
    t.canonical += ")";
    for (const auto& item: multiplicity)
        for (int j=2;j<=item.second;j++) t.automorphisms *= j;
    trees.push_back(t);
    return static_cast<int>(trees.size())-1;
}

struct Copy {
    int tree, side, offset;
    std::vector<int> positions;
};

struct Component {
    int tree, side; // 0: F side, 1: selected child, 2: H side.
    std::vector<int> positions;
};

struct Case {
    std::string name;
    bool energy;
    std::vector<int> f,h;
};

struct Audit {
    Case test;
    std::vector<Copy> copies;
    std::vector<int> owner, local, partner, block;
    std::vector<std::pair<int,int>> edges;
    int spin_count=0, special=-1;
    U64 admissible=0, empty=0, structural=0, disagreement=0;
    U64 non_eulerian_moment=0, no_free_edge_moment=0, non_tree_empty=0;
    U64 candidate_pairings=0;
    std::vector<int> first_bad;

    explicit Audit(Case c):test(std::move(c)) {
        auto append = [&](int tree, int side) {
            Copy cp{tree,side,spin_count,{}};
            int index = static_cast<int>(copies.size());
            const auto& parents = trees[tree].parent;
            for (int v=0;v<static_cast<int>(parents.size());v++) {
                cp.positions.push_back(spin_count);
                owner.push_back(index); local.push_back(v); spin_count++;
            }
            for (int v=0;v<static_cast<int>(parents.size());v++) {
                int endpoint = parents[v] < 0 ? (side==0 ? -1 : -2) : cp.offset+parents[v];
                edges.emplace_back(cp.offset+v,endpoint);
            }
            copies.push_back(cp);
        };
        for (int t:test.f) append(t,0);
        for (int t:test.h) append(t,1);
        if (test.energy) {
            special=spin_count++;
            owner.push_back(-1); local.push_back(-1);
            edges.emplace_back(-1,-2);
        }
        partner.assign(spin_count,-1);
        block.assign(spin_count,-1);
        if (spin_count%2) throw std::runtime_error("odd spin-count case");
    }

    bool allowed(int a,int b) const {
        if (owner[a]>=0 && owner[a]==owner[b]) return false;
        // An H copy is injective together with its root j, represented by
        // the explicit special spin. This restriction is independent of
        // where a surviving bridge partner would later have to lie.
        if (a==special || b==special) {
            int other = a==special ? b : a;
            if (owner[other]>=0 && copies[owner[other]].side==1) return false;
        }
        return true;
    }

    bool component_pairing(const std::vector<Component>& components) const {
        std::vector<int> which(spin_count,-1), index(spin_count,-1);
        for (int c=0;c<static_cast<int>(components.size());c++)
            for (int j=0;j<static_cast<int>(components[c].positions.size());j++) {
                int p=components[c].positions[j]; which[p]=c; index[p]=j;
            }
        for (int c=0;c<static_cast<int>(components.size());c++) {
            const auto& comp=components[c];
            int d=which[partner[comp.positions[0]]];
            if (d<0 || d==c) return false;
            const auto& other=components[d];
            if (trees[comp.tree].canonical!=trees[other.tree].canonical) return false;
            if (comp.side==0 || other.side==0) {
                if (comp.side!=0 || other.side!=0) return false;
            } else if (comp.side==1 && other.side==1) return false;
            const auto& par=trees[comp.tree].parent;
            const auto& other_par=trees[other.tree].parent;
            for (int j=0;j<static_cast<int>(comp.positions.size());j++) {
                int mate=partner[comp.positions[j]];
                if (which[mate]!=d) return false;
                int target=index[mate];
                if (par[j]<0) {
                    if (other_par[target]>=0) return false;
                } else {
                    int mate_parent=partner[comp.positions[par[j]]];
                    if (other_par[target]!=index[mate_parent]) return false;
                }
            }
        }
        return true;
    }

    bool is_structural() const {
        std::vector<Component> comp;
        if (!test.energy) {
            for (const auto& cp:copies) comp.push_back({cp.tree,0,cp.positions});
            return component_pairing(comp);
        }
        int top=partner[special];
        int selected=owner[top];
        if (selected<0 || copies[selected].side!=0 || local[top]!=0) return false;
        for (int c=0;c<static_cast<int>(copies.size());c++) {
            if (c==selected) continue;
            const auto& cp=copies[c];
            comp.push_back({cp.tree,cp.side==0 ? 0 : 2,cp.positions});
        }
        const Copy& selected_copy=copies[selected];
        int offset=selected_copy.offset+1;
        for (int child:trees[selected_copy.tree].children) {
            std::vector<int> positions;
            for (int j=0;j<static_cast<int>(trees[child].parent.size());j++)
                positions.push_back(offset+j);
            comp.push_back({child,1,positions});
            offset+=static_cast<int>(trees[child].parent.size());
        }
        return component_pairing(comp);
    }

    void inspect() {
        admissible++;
        int next=0;
        for (int p=0;p<spin_count;p++) if (p<partner[p])
            block[p]=block[partner[p]]=next++;
        const int root=next, vertices=next+1;
        int count[16][16] = {};
        auto label = [&](int p) {
            if (p==-1) return root;
            if (p==-2) return block[special];
            return block[p];
        };
        for (auto e:edges) {
            int a=label(e.first),b=label(e.second);
            if (a==b) throw std::runtime_error("injectivity missed a loop");
            if (a>b) std::swap(a,b);
            count[a][b]++;
        }
        bool zero=true, free_edge=false;
        int degree[16] = {}, distinct=0;
        for (int a=0;a<vertices;a++) for (int b=a+1;b<vertices;b++) {
            if (count[a][b]) distinct++;
            if (count[a][b]%2) {
                zero=false; degree[a]++; degree[b]++;
                if (a!=root && b!=root) free_edge=true;
            }
        }
        bool pattern=is_structural();
        empty+=zero; structural+=pattern;
        if (zero!=pattern) {
            disagreement++;
            if (first_bad.empty()) first_bad=partner;
        }
        if (!test.energy && !zero) {
            for (int d:degree) if (d%2) {non_eulerian_moment++; break;}
            if (!free_edge) no_free_edge_moment++;
        }
        if (zero) {
            bool bad=distinct!=vertices-1;
            for (int a=0;a<vertices;a++) for (int b=a+1;b<vertices;b++)
                if (count[a][b] && count[a][b]!=2) bad=true;
            if (bad) non_tree_empty++;
        }
    }

    void enumerate(std::uint32_t unused) {
        if (!unused) {inspect();return;}
        int a=__builtin_ctz(unused);
        std::uint32_t rest=unused & ~(1u<<a);
        for (std::uint32_t options=rest;options;options&=options-1) {
            int b=__builtin_ctz(options);
            candidate_pairings++;
            if (!allowed(a,b)) continue;
            partner[a]=b; partner[b]=a;
            enumerate(rest & ~(1u<<b));
        }
    }
};

U64 weighted_pairings(std::vector<std::pair<int,bool>> items) {
    if (items.empty()) return 1;
    U64 result=0;
    for (int j=1;j<static_cast<int>(items.size());j++) {
        if (items[0].first!=items[j].first || (items[0].second && items[j].second)) continue;
        std::vector<std::pair<int,bool>> remaining;
        for (int k=1;k<static_cast<int>(items.size());k++) if (k!=j) remaining.push_back(items[k]);
        result+=trees[items[0].first].automorphisms*weighted_pairings(remaining);
    }
    return result;
}

U64 prediction(const Case& c) {
    if (!c.energy) {
        std::vector<std::pair<int,bool>> items;
        for (int t:c.f) items.emplace_back(t,false);
        return weighted_pairings(items);
    }
    U64 answer=0;
    for (int selected=0;selected<static_cast<int>(c.f.size());selected++) {
        std::vector<std::pair<int,bool>> left,right;
        for (int j=0;j<static_cast<int>(c.f.size());j++) if (j!=selected) left.emplace_back(c.f[j],false);
        for (int t:trees[c.f[selected]].children) right.emplace_back(t,true);
        for (int t:c.h) right.emplace_back(t,false);
        answer+=weighted_pairings(left)*weighted_pairings(right);
    }
    return answer;
}

int main() {
    int e=add_tree("edge",{});
    int s=add_tree("star3",{e,e});
    int a5=add_tree("asymmetric5",{e,s});
    int a7=add_tree("asymmetric7",{e,a5});
    int b7=add_tree("balanced7",{s,s});
    std::vector<Case> cases={
        {"moment_edge4",false,{e,e,e,e},{}},
        {"moment_star3_4",false,{s,s,s,s},{}},
        {"moment_asymmetric5_2",false,{a5,a5},{}},
        {"moment_asymmetric7_2",false,{a7,a7},{}},
        {"moment_balanced7_2",false,{b7,b7},{}},
        {"moment_mixed_1_3_5_7",false,{e,s,a5,a7},{}},
        {"moment_star3_2_asymmetric5_2",false,{s,s,a5,a5},{}},
        {"moment_asymmetric7_2_edge2",false,{a7,a7,e,e},{}},
        {"moment_balanced7_2_edge2",false,{b7,b7,e,e},{}},
        {"moment_edge4_star3_2",false,{e,e,e,e,s,s},{}},
        {"moment_asymmetric5_4",false,{a5,a5,a5,a5},{}},
        {"moment_balanced7_2_star3_2",false,{b7,b7,s,s},{}},
        {"moment_mixed_7_7_5_1",false,{a7,b7,a5,e},{}},
        {"energy_edge_empty",true,{e},{}},
        {"energy_star3_edge2",true,{s},{e,e}},
        {"energy_asymmetric5_edge_star3",true,{a5},{e,s}},
        {"energy_asymmetric7_edge_asymmetric5",true,{a7},{e,a5}},
        {"energy_balanced7_star3_2",true,{b7},{s,s}},
        {"energy_asymmetric5_edge3_star3",true,{a5},{e,e,e,s}},
        {"energy_edge_star3_2_empty",true,{e,s,s},{}},
        {"energy_star3_3_edge2",true,{s,s,s},{e,e}},
        {"energy_star3_asymmetric5_2_edge2",true,{s,a5,a5},{e,e}},
        {"energy_balanced7_star3_2_edge2",true,{b7},{s,s,e,e}},
        {"energy_asymmetric7_asymmetric5_edge3",true,{a7},{a5,e,e,e}},
        {"energy_wrong_children_asymmetric7",true,{a7},{s,s}},
        {"energy_wrong_children_balanced7",true,{b7},{e,a5}},
        {"energy_even_F",true,{e,e},{e}},
        {"energy_mixed_zero",true,{a5,s,e},{s,e}},
        {"energy_asymmetric5_3_edge_star3",true,{a5,a5,a5},{e,s}},
    };
    bool pass=true;
    U64 total=0;
    std::cout<<"{\n  \"method\": \"exhaustive_exact_leading_pair_partitions\",\n  \"trees\": [\n";
    for (int j=0;j<static_cast<int>(trees.size());j++) {
        const auto& t=trees[j];
        std::cout<<"    {\"name\":\""<<t.name<<"\",\"nonroot_vertices\":"<<t.parent.size()
                 <<",\"automorphisms\":"<<t.automorphisms<<",\"canonical\":\""<<t.canonical<<"\"}"
                 <<(j+1==static_cast<int>(trees.size()) ? "\n" : ",\n");
    }
    std::cout<<"  ],\n  \"cases\": [\n";
    for (int j=0;j<static_cast<int>(cases.size());j++) {
        Audit audit(cases[j]);
        if (audit.spin_count>20) throw std::runtime_error("case too large");
        audit.enumerate((1u<<audit.spin_count)-1);
        U64 expected=prediction(cases[j]);
        bool good=audit.disagreement==0 && audit.non_eulerian_moment==0
                  && audit.no_free_edge_moment==0 && audit.non_tree_empty==0
                  && audit.empty==expected;
        pass=pass && good; total+=audit.admissible;
        std::cout<<"    {\"name\":\""<<cases[j].name<<"\",\"spin_slots\":"<<audit.spin_count
                 <<",\"admissible_partitions\":"<<audit.admissible
                 <<",\"zero_parity_partitions\":"<<audit.empty
                 <<",\"structural_partitions\":"<<audit.structural
                 <<",\"independent_Wick_Hermite_prediction\":"<<expected
                 <<",\"classification_disagreements\":"<<audit.disagreement
                 <<",\"non_eulerian_moment_parities\":"<<audit.non_eulerian_moment
                 <<",\"moment_parities_without_free_edge\":"<<audit.no_free_edge_moment
                 <<",\"empty_parities_not_double_edge_trees\":"<<audit.non_tree_empty
                 <<",\"passed\":"<<(good ? "true" : "false")<<"}"
                 <<(j+1==static_cast<int>(cases.size()) ? "\n" : ",\n");
        if (!audit.first_bad.empty()) {
            std::cerr<<"COUNTEREXAMPLE "<<cases[j].name<<":";
            for (int p:audit.first_bad) std::cerr<<" "<<p;
            std::cerr<<"\n";
        }
    }
    std::cout<<"  ],\n  \"total_admissible_partitions\": "<<total
             <<",\n  \"verified\": "<<(pass ? "true" : "false")<<"\n}\n";
    return pass ? 0 : 1;
}
