test_that("Spherical harmonics", {
    local_edition(3)
    
    m <- inla.mesh.create(globe = 1)
    sph0 <- inla.mesh.basis(m, type = "sph.harm", rot.inv = TRUE)
    
    expect_equal(dim(sph0), c(12, 4))
    expect_equal(sph0[1, ], c(1, 1.732050808, 2.2360679775, 2.645751311))
    
    sph <- inla.mesh.basis(m, type = "sph.harm", rot.inv = FALSE)
    
    expect_equal(dim(sph), c(12, 16))
    expect_equal(sph[1, seq_len(4)], c(1, 0 , 1.732050808, 0))
    expect_equal(sph[1, c(1, 3, 7, 13)],
                 c(1, 1.732050808, 2.2360679775, 2.645751311))
})
