#!/bin/bash

runNumberList=( 071728 071729 071734 071735 071736 071737 071739 071740 071742 071743 071747 071748 071751 071752 071753 071754 071755 071757 071758 071759 )
for i in "${runNumberList[@]}"
do
   root -b -q 'TAGMTWExtractor.C("'$i'")' &
done

# Not running yet
# 071728 071729 071734 071735 071736 071737 071739 071740 071742 071743 071747 071748 071751 071752 071753 071754 071755 071757 071758 071759
# 071760 071761 071762 071763 071765 071766 071767 071768 071769 071770 071772 071773 071774 071775 071776 071777 071778 071780 071781 071782 071783 071784
# 071785 071786 071789 071790 071791 071792 071793 071794 071795 071800 071806 071807 071808 071809
# 071810 071813 071814 071815 071817 071819 071820 071821 071822 071823 071824 071825 071826 071827
# 071830 071831 071832 071833 071834 071835 071836 071837 071838 071839 071840 071841 071842 071846


#Completed
# 071847 071848 071849 071850 071851 071852 071853 071854 071855 071863 071864 071865 071867 071868 071869 071870
