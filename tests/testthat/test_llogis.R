context("Log-logistic distribution")

test_that("llogis",{
    x <- c(0.1, 0.2, 0.7)
    if (require("eha"))
        expect_equal(dllogis(x, shape=0.1, scale=0.2), eha::dllogis(x, shape=0.1, scale=0.2))
    expect_equal(qllogis(x, shape=0.1, scale=0.2), qgeneric(pllogis, p=x, shape=0.1, scale=0.2))
    expect_equal(x, pllogis(qllogis(x, shape=0.1, scale=0.2), shape=0.1, scale=0.2))
    expect_equal(x, 1 - exp(-Hllogis(qllogis(x, shape=0.1, scale=0.2), shape=0.1, scale=0.2)))
    expect_equal(hllogis(x, shape=0.1, scale=0.2),  dllogis(x, shape=0.1, scale=0.2) / (1 - pllogis(x, shape=0.1, scale=0.2)))
    q <- numeric(3); for (i in 1:3) q[i] <- qllogis(x[i], shape=0.0001, scale=i/5)
    expect_equal(q, qllogis(x, shape=0.0001, scale=1:3/5))
    x <- c(0.5, 1.06, 4.7)
    expect_equal(x, qllogis(pllogis(x, shape=0.1, scale=0.2), shape=0.1, scale=0.2))
    if (interactive()) {
        rl <- rllogis(10000, shape=1.5, scale=1.2)
        plot(density(rl[rl<100]), xlim=c(0,10))
        x <- seq(0, 10, by=0.001)
        lines(x, dllogis(x, shape=1.5, scale=1.2), lty=2)
    }
    expect_equal(mean_llogis(shape=0.1, scale=0.2), NaN)
    expect_equal(var.llogis(shape=1.1, scale=0.2), NaN)
    mean_llogis(shape=1.1, scale=0.2)
    var.llogis(shape=2.1, scale=0.2)

    x <- c(0.1,0.4,0.6,0.9)
    expect_equal(pllogis(qllogis(p=x, lower.tail=FALSE), lower.tail=FALSE), x)
    expect_equal(qlogis(p=x, lower.tail=TRUE), rev(qlogis(p=x, lower.tail=FALSE)))
    
})
