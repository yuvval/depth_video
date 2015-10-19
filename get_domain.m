function domain_alias = get_domain()


[~, domain] = system('dnsdomainname');

if is_in_str('lnx.biu.ac.il', domain)
    domain_alias = 'cortex';
elseif is_in_str('csail.mit.edu', domain)
    domain_alias = 'csail';
end

